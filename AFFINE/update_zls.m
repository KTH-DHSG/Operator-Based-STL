function update_zls(plhandle,zls_state_data,zls_time_data,treach,tend,Color)

    plhandle.XData = [tend-treach,tend-treach,zls_time_data{1},nan,tend-treach,zls_time_data{2}] + treach;

    plhandle.YData = [zls_state_data{2}(1),zls_state_data{1}(1),zls_state_data{1},nan,zls_state_data{2}(1),zls_state_data{2}];

    plhandle.Color = Color;
end