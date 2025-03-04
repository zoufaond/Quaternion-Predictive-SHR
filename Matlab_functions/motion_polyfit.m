function [time_second, data_second]= motion_polyfit(motion_full,first_interp_num,second_interp_num,plot_polyfits)

time = motion_full(:,1)-motion_full(1,1);
time_first = linspace(0,time(end,1),first_interp_num)';
time_second = linspace(0,time(end,1),second_interp_num)';
mot_eul = deg2rad(motion_full(:,8:end)); %ignore time and rotation/translation of thorax
data_first = zeros(first_interp_num,11);
data_second = zeros(second_interp_num,11);

for jnt = 1:length(mot_eul(1,:))
    data_first(:,jnt) = interp1(time,mot_eul(:,jnt),time_first,'spline');
    data_second(:,jnt) = interp1(time_first,data_first(:,jnt),time_second,'spline');
    if plot_polyfits == 1
        figure
        plot(time,mot_eul(:,jnt),'o')
        hold on
        plot(time_first,data_first(:,jnt))
        hold on 
        plot(time_second,data_second(:,jnt))
        legend('IK','First interp','Second_interp')
    end
end

end