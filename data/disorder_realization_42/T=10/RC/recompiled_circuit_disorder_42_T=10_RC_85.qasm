OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5126098) q[0];
sx q[0];
rz(-1.0273758) q[0];
sx q[0];
rz(0.36261121) q[0];
rz(0.91712046) q[1];
sx q[1];
rz(2.6511104) q[1];
sx q[1];
rz(9.766415) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5255614) q[0];
sx q[0];
rz(-1.6934868) q[0];
sx q[0];
rz(0.37567715) q[0];
rz(-pi) q[1];
rz(0.59092893) q[2];
sx q[2];
rz(-1.8701631) q[2];
sx q[2];
rz(0.17609827) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.29400723) q[1];
sx q[1];
rz(-2.5348186) q[1];
sx q[1];
rz(-2.1104382) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2239286) q[3];
sx q[3];
rz(-1.9609946) q[3];
sx q[3];
rz(1.8412561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.59387702) q[2];
sx q[2];
rz(-1.1527858) q[2];
sx q[2];
rz(0.95648742) q[2];
rz(2.9603413) q[3];
sx q[3];
rz(-2.5528208) q[3];
sx q[3];
rz(1.6199002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0797121) q[0];
sx q[0];
rz(-0.065608874) q[0];
sx q[0];
rz(1.99615) q[0];
rz(2.0727797) q[1];
sx q[1];
rz(-0.83199465) q[1];
sx q[1];
rz(2.4157445) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94978588) q[0];
sx q[0];
rz(-1.0254142) q[0];
sx q[0];
rz(-1.9907065) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5475694) q[2];
sx q[2];
rz(-1.1813287) q[2];
sx q[2];
rz(2.4246852) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6594071) q[1];
sx q[1];
rz(-0.30989753) q[1];
sx q[1];
rz(1.6531497) q[1];
rz(-pi) q[2];
rz(-1.3594567) q[3];
sx q[3];
rz(-2.1826715) q[3];
sx q[3];
rz(-0.26821995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6590865) q[2];
sx q[2];
rz(-1.7616452) q[2];
sx q[2];
rz(-2.144311) q[2];
rz(1.4128134) q[3];
sx q[3];
rz(-1.0935067) q[3];
sx q[3];
rz(-2.2028082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41855758) q[0];
sx q[0];
rz(-0.96590531) q[0];
sx q[0];
rz(-1.1285271) q[0];
rz(1.8380802) q[1];
sx q[1];
rz(-2.4120657) q[1];
sx q[1];
rz(-0.9054786) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.758924) q[0];
sx q[0];
rz(-1.3805458) q[0];
sx q[0];
rz(2.8959136) q[0];
x q[1];
rz(0.75240527) q[2];
sx q[2];
rz(-1.89626) q[2];
sx q[2];
rz(2.4706555) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14721522) q[1];
sx q[1];
rz(-2.4265687) q[1];
sx q[1];
rz(1.6836402) q[1];
rz(2.2842992) q[3];
sx q[3];
rz(-1.6033353) q[3];
sx q[3];
rz(-1.1710222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3811615) q[2];
sx q[2];
rz(-0.16813777) q[2];
sx q[2];
rz(-2.1543489) q[2];
rz(3.0958214) q[3];
sx q[3];
rz(-1.8493303) q[3];
sx q[3];
rz(-0.18019095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0320597) q[0];
sx q[0];
rz(-0.68234545) q[0];
sx q[0];
rz(3.0155244) q[0];
rz(2.9571422) q[1];
sx q[1];
rz(-2.0729063) q[1];
sx q[1];
rz(1.9741845) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084598736) q[0];
sx q[0];
rz(-1.3697249) q[0];
sx q[0];
rz(0.76954557) q[0];
x q[1];
rz(-1.3912958) q[2];
sx q[2];
rz(-1.0738392) q[2];
sx q[2];
rz(2.3617982) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0988401) q[1];
sx q[1];
rz(-2.4272857) q[1];
sx q[1];
rz(0.40978281) q[1];
rz(-0.79951841) q[3];
sx q[3];
rz(-2.3164253) q[3];
sx q[3];
rz(-1.0471917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2771153) q[2];
sx q[2];
rz(-2.7272759) q[2];
sx q[2];
rz(-2.5100822) q[2];
rz(-0.19550066) q[3];
sx q[3];
rz(-1.86444) q[3];
sx q[3];
rz(0.41803944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1642078) q[0];
sx q[0];
rz(-3.003484) q[0];
sx q[0];
rz(-0.50267977) q[0];
rz(-2.0282822) q[1];
sx q[1];
rz(-2.138442) q[1];
sx q[1];
rz(-0.46708333) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57468092) q[0];
sx q[0];
rz(-1.7705288) q[0];
sx q[0];
rz(-1.9499705) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3181114) q[2];
sx q[2];
rz(-0.99950302) q[2];
sx q[2];
rz(-0.90734446) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0568697) q[1];
sx q[1];
rz(-1.1493756) q[1];
sx q[1];
rz(0.23878581) q[1];
rz(1.763837) q[3];
sx q[3];
rz(-0.69460624) q[3];
sx q[3];
rz(-0.65929268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.46665329) q[2];
sx q[2];
rz(-1.8920687) q[2];
sx q[2];
rz(-0.49218991) q[2];
rz(0.28218937) q[3];
sx q[3];
rz(-1.4331093) q[3];
sx q[3];
rz(1.5757489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.870134) q[0];
sx q[0];
rz(-2.6014355) q[0];
sx q[0];
rz(3.0555994) q[0];
rz(-1.7135235) q[1];
sx q[1];
rz(-1.0079039) q[1];
sx q[1];
rz(-1.6451947) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2841227) q[0];
sx q[0];
rz(-0.5623445) q[0];
sx q[0];
rz(1.9447295) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9005269) q[2];
sx q[2];
rz(-2.0965555) q[2];
sx q[2];
rz(-0.17503967) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.306957) q[1];
sx q[1];
rz(-0.35528696) q[1];
sx q[1];
rz(-1.5026232) q[1];
x q[2];
rz(-1.5835285) q[3];
sx q[3];
rz(-2.2836962) q[3];
sx q[3];
rz(1.7604148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4795586) q[2];
sx q[2];
rz(-0.61926121) q[2];
sx q[2];
rz(-0.39624828) q[2];
rz(2.5475492) q[3];
sx q[3];
rz(-2.0500573) q[3];
sx q[3];
rz(-1.5007639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632161) q[0];
sx q[0];
rz(-2.5686869) q[0];
sx q[0];
rz(1.0923882) q[0];
rz(1.3573525) q[1];
sx q[1];
rz(-1.4211632) q[1];
sx q[1];
rz(3.1299652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25225885) q[0];
sx q[0];
rz(-0.47124915) q[0];
sx q[0];
rz(-2.528119) q[0];
rz(-pi) q[1];
rz(1.5426703) q[2];
sx q[2];
rz(-1.5962432) q[2];
sx q[2];
rz(-0.59288073) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1639451) q[1];
sx q[1];
rz(-1.5944591) q[1];
sx q[1];
rz(1.0246828) q[1];
x q[2];
rz(-1.8947951) q[3];
sx q[3];
rz(-0.9231336) q[3];
sx q[3];
rz(-0.049829359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.45550436) q[2];
sx q[2];
rz(-0.65767613) q[2];
sx q[2];
rz(-0.28798506) q[2];
rz(1.1503495) q[3];
sx q[3];
rz(-1.8201613) q[3];
sx q[3];
rz(-0.59818017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8312254) q[0];
sx q[0];
rz(-2.6103525) q[0];
sx q[0];
rz(1.2121375) q[0];
rz(0.37777004) q[1];
sx q[1];
rz(-1.9394082) q[1];
sx q[1];
rz(-1.6479962) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15570116) q[0];
sx q[0];
rz(-2.771286) q[0];
sx q[0];
rz(-0.62472384) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0817501) q[2];
sx q[2];
rz(-1.4388196) q[2];
sx q[2];
rz(3.0795385) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6495325) q[1];
sx q[1];
rz(-2.3146555) q[1];
sx q[1];
rz(2.5745113) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0231608) q[3];
sx q[3];
rz(-1.0716972) q[3];
sx q[3];
rz(-2.0407608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7887855) q[2];
sx q[2];
rz(-2.0264503) q[2];
sx q[2];
rz(1.8898233) q[2];
rz(-2.2552323) q[3];
sx q[3];
rz(-0.76615196) q[3];
sx q[3];
rz(-0.10929508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66824526) q[0];
sx q[0];
rz(-1.0412403) q[0];
sx q[0];
rz(-2.9633203) q[0];
rz(-0.072862236) q[1];
sx q[1];
rz(-2.5477414) q[1];
sx q[1];
rz(-1.1791621) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0807063) q[0];
sx q[0];
rz(-1.0651677) q[0];
sx q[0];
rz(1.9067184) q[0];
x q[1];
rz(0.945325) q[2];
sx q[2];
rz(-1.1497467) q[2];
sx q[2];
rz(-0.23907121) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0563911) q[1];
sx q[1];
rz(-1.8279652) q[1];
sx q[1];
rz(-0.53253865) q[1];
rz(-pi) q[2];
rz(-0.9411373) q[3];
sx q[3];
rz(-0.59560532) q[3];
sx q[3];
rz(-2.5411118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4420085) q[2];
sx q[2];
rz(-2.1506491) q[2];
sx q[2];
rz(1.0317831) q[2];
rz(-1.7992841) q[3];
sx q[3];
rz(-1.4167901) q[3];
sx q[3];
rz(0.22527307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28829065) q[0];
sx q[0];
rz(-1.0422215) q[0];
sx q[0];
rz(3.0798966) q[0];
rz(1.8348947) q[1];
sx q[1];
rz(-1.6625762) q[1];
sx q[1];
rz(1.0888938) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0058044) q[0];
sx q[0];
rz(-3.0514768) q[0];
sx q[0];
rz(0.80149217) q[0];
rz(-pi) q[1];
rz(2.9577272) q[2];
sx q[2];
rz(-2.2176748) q[2];
sx q[2];
rz(1.136214) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1354462) q[1];
sx q[1];
rz(-1.0137614) q[1];
sx q[1];
rz(-2.1250493) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0652566) q[3];
sx q[3];
rz(-2.7571602) q[3];
sx q[3];
rz(-1.2151375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5237727) q[2];
sx q[2];
rz(-0.68796316) q[2];
sx q[2];
rz(3.0715959) q[2];
rz(2.2648515) q[3];
sx q[3];
rz(-1.3249967) q[3];
sx q[3];
rz(2.783412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6488279) q[0];
sx q[0];
rz(-1.5179101) q[0];
sx q[0];
rz(-2.5483325) q[0];
rz(1.5169253) q[1];
sx q[1];
rz(-2.0943874) q[1];
sx q[1];
rz(2.692683) q[1];
rz(-1.0629874) q[2];
sx q[2];
rz(-0.72401902) q[2];
sx q[2];
rz(0.25821092) q[2];
rz(-0.041257507) q[3];
sx q[3];
rz(-1.9058803) q[3];
sx q[3];
rz(-2.1432721) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
