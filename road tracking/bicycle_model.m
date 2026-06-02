function xnext = bicycle_model(x, w)
    L = 0.137;
    dt = 0.5;
    
    px = x(1);
    py = x(2);
    theta = x(3);
    v = x(4);
    
    a_noise = w(3);
    delta_noise = w(4);
    
    px_next = px + v*cos(theta)*dt+w(1)*dt;
    py_next = py + v*sin(theta)*dt+w(2)*dt;
    theta_next = theta + delta_noise*dt;
    v_next = v + a_noise*dt;
    
    xnext = [px_next; py_next; theta_next; v_next];
end
