configs = [noloadlin,noloadnonlin,zerofieldlin,zerofieldnonlin,...
    loadlin,loadnonlin];
names = ["brad_noload_lin","brad_noload_nonlin","brad_zerofield_lin",...
    "brad_zerofield_nonlin","brad_load_lin","brad_load_nonlin"];

for i = 1:size(configs,2)
    figure
    plot(angle,configs(:,i))
    xlabel("Angle [rad]")
    ylabel("Radial flux [Wb]")
    xlim([0,angle(end)])
    print(names(i),'-dpng','-r500')
end