export nondivergent_vector_field_2d, make_stream_function

function stream_function(xX, yY, grid, dx, dy, omega, xc, yc)
    x = xX * grid[1] * dx
    y = yY * grid[2] * dy
    return 0.5 * omega * ((x .- xc).^2 .+ (y .- yc).^2)
end

function make_stream_function(grid, dx, dy, omega, xc, yc)
    sf = (x, y) -> stream_function(x, y, grid, dx, dy, omega, xc, yc)
    return sf
end

function face_coords(grid, offset_x=0.0, offset_y=0.0)
    nx, ny = grid
    if offset_x > 0
        ny += 1
    elseif offset_y > 0
        nx += 1
    end
    x = range(offset_x, stop=nx - 1 + offset_x, length=nx) ./ nx
    y = range(offset_y, stop=ny - 1 + offset_y, length=ny) ./ ny
    xX = repeat(x, 1, ny)
    yY = repeat(collect(y)', nx, 1)
    return xX, yY
end

function nondivergent_vector_field_2d(dx,dy,grid, dt, sf)
    dxX, dyY = 1 ./ grid

    xX, yY = face_coords(grid, 0.0, 0.5)
    u_velocity = -1 .* (sf(xX, yY .+ dyY / 2) .- sf(xX, yY .- dyY / 2)) ./ dy

    xX, yY = face_coords(grid, 0.5, 0.0)
    v_velocity = (sf(xX .+ dxX / 2, yY) .- sf(xX .- dxX / 2, yY)) ./ dx

    GCx = u_velocity * dt / dx
    GCy = v_velocity * dt / dy

    @assert all(abs.(GCx) .< 1.0) "CFL violation in x"
    @assert all(abs.(GCy) .< 1.0) "CFL violation in y"

    return u_velocity, v_velocity,GCx,GCy
end