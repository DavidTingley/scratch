
clear
d  = dir('*201*');
wheel_count = 1;  % counters for trial numbers of each type
central_count = 1;
linear_count = 1;
 
%% compile data

for ii=1:length(d)
   cd(d(ii).name)
   load([d(ii).name '.behavior.mat'],'behavior')
   if strcmp(behavior.units,'pixels')  % convert to cm for all tracking types...
       multiplier = 6.5;
   elseif strcmp(behavior.units,'m')
       multiplier = 10000;
   elseif strcmp(behavior.units,'mm')
       multiplier = 10;
   else
       error('could not find correct multiplier')
   end

   if strcmp(behavior.description,'wheel alternation')
       for trial = 1:length(behavior.events.trials)
           wheel_vel(wheel_count,:) = medfilt1(makeLength(sqrt((abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))).^2),200),20) .* multiplier;
           wheel_time(wheel_count) = behavior.events.trials{trial}.timestamps(end) - behavior.events.trials{trial}.timestamps(1);
           wheel_units{wheel_count} = behavior.units;
           if strcmp(behavior.units,'pixels')
               wheel_dist(wheel_count) = sum(abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))) * multiplier * 3.9;
           else
               wheel_dist(wheel_count) = sum(abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))) * multiplier;
           end
           if sum(wheel_vel(wheel_count,:))<100
              error() 
           end
           wheel_count = wheel_count+1;
       end
   elseif strcmp(behavior.description,'central alternation')
        for trial = 1:length(behavior.events.trials)
            central_vel(central_count,:) = medfilt1(makeLength(sqrt((abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))).^2),200),20) .* multiplier;
            central_time(central_count) = behavior.events.trials{trial}.timestamps(end) - behavior.events.trials{trial}.timestamps(1);
            central_units{central_count} = behavior.units;
            central_dist(central_count) = sum(abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))) * multiplier;
            if sum(central_vel(central_count,:))<100
            error() 
            end
                       
                       central_count = central_count+1;
        end
   elseif strcmp(behavior.description,'linear')
        for trial = 1:length(behavior.events.trials)
            linear_vel(linear_count,:) = medfilt1(makeLength(sqrt((abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))).^2),200),20) .* multiplier;
            linear_time(linear_count) = behavior.events.trials{trial}.timestamps(end) - behavior.events.trials{trial}.timestamps(1);
            linear_units{linear_count} = behavior.units;
            linear_dist(linear_count) = sum(abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))) * multiplier;
            if sum(linear_vel(linear_count,:))<100
            error() 
            end
                        linear_count = linear_count+1;
        end
   elseif strcmp(behavior.description,'both alternation')
       for condition = 1:length(unique(behavior.events.trialConditions))
          % determine condition type...
          f = find(behavior.events.trialConditions==condition);
          for trial = 1:length(f)
            dist(f) = sum(abs(diff(behavior.events.trials{f(trial)}.x)) + abs(diff(behavior.events.trials{f(trial)}.y)));
          end
          [a b] = min([35278 29743] - mean(dist));
          % assign data
          if b == 1
               for t = 1:length(f)
                   trial = f(t);
                   wheel_vel(wheel_count,:) = medfilt1(makeLength(sqrt((abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))).^2),200),20) .* multiplier;
                   wheel_time(wheel_count) = behavior.events.trials{trial}.timestamps(end) - behavior.events.trials{trial}.timestamps(1);
                   wheel_units{wheel_count} = behavior.units;
                   if strcmp(behavior.units,'pixels')
                       wheel_dist(wheel_count) = sum(abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))) * multiplier * 3.9;
                   else
                       wheel_dist(wheel_count) = sum(abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))) * multiplier;
                   end
                    if sum(wheel_vel(wheel_count,:))<100
                    error() 
                    end
                   wheel_count = wheel_count+1;
               end
          elseif b == 2
            for t = 1:length(f)
                trial = f(t);
                central_vel(central_count,:) = medfilt1(makeLength(sqrt((abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))).^2),200),20) .* multiplier;
                central_time(central_count) = behavior.events.trials{trial}.timestamps(end) - behavior.events.trials{trial}.timestamps(1);
                central_units{central_count} = behavior.units;
                central_dist(central_count) = sum(abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))) * multiplier;
                if sum(central_vel(central_count,:))<100
                error() 
                end
                            central_count = central_count+1;
            end
          end
       end
   elseif strcmp(behavior.description,'linear/jump')
       for trial = 1:length(behavior.events.trials)
          if strcmp(behavior.events.conditionType{behavior.events.trialConditions(trial)},'linear')
            linear_vel(linear_count,:) = medfilt1(makeLength(sqrt((abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))).^2),200),20) .* multiplier;
            linear_time(linear_count) = behavior.events.trials{trial}.timestamps(end) - behavior.events.trials{trial}.timestamps(1);
            linear_units{linear_count} = behavior.units;
            linear_dist(linear_count) = sum(abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))) * multiplier;
            linear_count = linear_count+1;
          elseif strcmp(behavior.events.conditionType{behavior.events.trialConditions(trial)},'central')
            central_vel(central_count,:) = medfilt1(makeLength(sqrt((abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))).^2),200),20) .* multiplier;
            central_time(central_count) = behavior.events.trials{trial}.timestamps(end) - behavior.events.trials{trial}.timestamps(1);
            central_units{central_count} = behavior.units;
            central_dist(central_count) = sum(abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))) * multiplier;
            central_count = central_count+1;
          elseif strcmp(behavior.events.conditionType{behavior.events.trialConditions(trial)},'wheel')
            wheel_vel(wheel_count,:) = medfilt1(makeLength(sqrt((abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))).^2),200),20) .* multiplier;
            wheel_time(wheel_count) = behavior.events.trials{trial}.timestamps(end) - behavior.events.trials{trial}.timestamps(1);
            wheel_units{wheel_count} = behavior.units;
            if strcmp(behavior.units,'pixels')
                wheel_dist(wheel_count) = sum(abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))) * multiplier * 3.9;
            else
                wheel_dist(wheel_count) = sum(abs(diff(behavior.events.trials{trial}.x)) + abs(diff(behavior.events.trials{trial}.y))) * multiplier;
            end
            wheel_count = wheel_count+1;
          else
             error 
          end
       end
   else
       error('couldnt find behavior')
   end
        cd /home/david/datasets/lsDataset/
end


%% start plotting


subplot(4,4,1)
boundedline(0:.01601:3.2,nanmean(wheel_vel),nanstd(wheel_vel))
axis([0 3.2 0 120])
ylabel('cm/s')
xlabel('position (meters)')
title('circle wheel alternation')
subplot(4,4,2)
histogram(wheel_vel(:),0:5:150)
xlabel('cm/s')
title('all velocities')
subplot(4,4,3)
histogram(mean(wheel_vel'),0:5:150)
xlabel('cm/s')
title('mean velocities')
subplot(4,4,4)
histogram(wheel_time(:),1:.2:10)
title('trial duration')
xlabel('time (seconds)')

subplot(4,4,5)
boundedline(0:.01301:2.6,nanmean(central_vel),nanstd(central_vel))
axis([0 2.6 0 120])
ylabel('cm/s')
xlabel('position (meters)')
title('central stem alternation')
subplot(4,4,6)
histogram(central_vel(:),0:5:150)
title('all velocities')
xlabel('cm/s')
subplot(4,4,7)
histogram(mean(central_vel'),0:5:150)
title('mean velocities')
xlabel('cm/s')
subplot(4,4,8)
histogram(central_time(:),1:.2:10)
title('trial duration')
xlabel('time (seconds)')

subplot(4,4,9)
boundedline(0:.01001:2,nanmean(central_vel),nanstd(central_vel))
axis([0 2 0 120])
ylabel('cm/s')
xlabel('position (meters)')
title('linear track')
subplot(4,4,10)
histogram(linear_vel(:),0:5:150)
title('all velocities')
xlabel('cm/s')
subplot(4,4,11)
histogram(mean(linear_vel'),0:5:150)
title('mean velocities')
xlabel('cm/s')
subplot(4,4,12)
histogram(linear_time(:),1:.2:10)
title('trial duration')
xlabel('time (seconds)')








