OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.81340462) q[0];
sx q[0];
rz(8.8153664) q[0];
sx q[0];
rz(9.463268) q[0];
rz(2.6961532) q[1];
sx q[1];
rz(-0.69094509) q[1];
sx q[1];
rz(3.0145187) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2949383) q[0];
sx q[0];
rz(-1.9401258) q[0];
sx q[0];
rz(-2.3890004) q[0];
x q[1];
rz(1.5655925) q[2];
sx q[2];
rz(-1.570911) q[2];
sx q[2];
rz(-0.031129908) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5161612) q[1];
sx q[1];
rz(-1.2417842) q[1];
sx q[1];
rz(2.8686868) q[1];
rz(-pi) q[2];
rz(2.1640763) q[3];
sx q[3];
rz(-1.1597871) q[3];
sx q[3];
rz(1.051595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8700063) q[2];
sx q[2];
rz(-0.25124696) q[2];
sx q[2];
rz(1.1146438) q[2];
rz(-2.8031269) q[3];
sx q[3];
rz(-1.7287799) q[3];
sx q[3];
rz(0.72845212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.912643) q[0];
sx q[0];
rz(-2.8680608) q[0];
sx q[0];
rz(-1.9563414) q[0];
rz(1.7456906) q[1];
sx q[1];
rz(-1.4811367) q[1];
sx q[1];
rz(0.11298583) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4596371) q[0];
sx q[0];
rz(-1.5636866) q[0];
sx q[0];
rz(-0.05592859) q[0];
x q[1];
rz(0.24893181) q[2];
sx q[2];
rz(-2.0023731) q[2];
sx q[2];
rz(-0.54655036) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7020018) q[1];
sx q[1];
rz(-1.6338305) q[1];
sx q[1];
rz(-1.3443771) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91199888) q[3];
sx q[3];
rz(-2.5781401) q[3];
sx q[3];
rz(-0.4215301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9129755) q[2];
sx q[2];
rz(-1.7824219) q[2];
sx q[2];
rz(1.011298) q[2];
rz(-0.060981123) q[3];
sx q[3];
rz(-1.1119548) q[3];
sx q[3];
rz(2.8673867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36534742) q[0];
sx q[0];
rz(-1.341935) q[0];
sx q[0];
rz(1.1439398) q[0];
rz(-1.2155608) q[1];
sx q[1];
rz(-2.2823915) q[1];
sx q[1];
rz(2.5124195) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7157756) q[0];
sx q[0];
rz(-1.5377147) q[0];
sx q[0];
rz(2.1685409) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10444) q[2];
sx q[2];
rz(-2.0757489) q[2];
sx q[2];
rz(-1.4940408) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0720318) q[1];
sx q[1];
rz(-1.6545418) q[1];
sx q[1];
rz(-2.778758) q[1];
rz(1.6225624) q[3];
sx q[3];
rz(-1.8300149) q[3];
sx q[3];
rz(1.0362877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.7367785) q[2];
sx q[2];
rz(-1.7958769) q[2];
sx q[2];
rz(0.2573615) q[2];
rz(1.2579873) q[3];
sx q[3];
rz(-1.6519929) q[3];
sx q[3];
rz(-2.6902698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8266206) q[0];
sx q[0];
rz(-0.74084145) q[0];
sx q[0];
rz(-1.524087) q[0];
rz(-1.357088) q[1];
sx q[1];
rz(-2.4391104) q[1];
sx q[1];
rz(-1.9125028) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8600677) q[0];
sx q[0];
rz(-1.901214) q[0];
sx q[0];
rz(-1.9343253) q[0];
rz(2.4357834) q[2];
sx q[2];
rz(-2.2459111) q[2];
sx q[2];
rz(1.0256888) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7732765) q[1];
sx q[1];
rz(-1.6267836) q[1];
sx q[1];
rz(0.30572388) q[1];
rz(-pi) q[2];
rz(0.96276729) q[3];
sx q[3];
rz(-1.9809147) q[3];
sx q[3];
rz(-1.6671163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.412753) q[2];
sx q[2];
rz(-1.7985901) q[2];
sx q[2];
rz(3.0573678) q[2];
rz(-1.6526875) q[3];
sx q[3];
rz(-3.0907349) q[3];
sx q[3];
rz(2.1648255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65234891) q[0];
sx q[0];
rz(-0.69559613) q[0];
sx q[0];
rz(0.034164567) q[0];
rz(1.4802406) q[1];
sx q[1];
rz(-2.7378597) q[1];
sx q[1];
rz(1.2108948) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04003018) q[0];
sx q[0];
rz(-1.2759325) q[0];
sx q[0];
rz(1.9676137) q[0];
x q[1];
rz(-1.239085) q[2];
sx q[2];
rz(-1.0704652) q[2];
sx q[2];
rz(-2.1608888) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3386061) q[1];
sx q[1];
rz(-0.44281964) q[1];
sx q[1];
rz(-1.3989596) q[1];
rz(-1.7922395) q[3];
sx q[3];
rz(-1.4998271) q[3];
sx q[3];
rz(-2.8896795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.5517129) q[2];
sx q[2];
rz(-0.61673841) q[2];
sx q[2];
rz(-1.7390772) q[2];
rz(-1.2276522) q[3];
sx q[3];
rz(-0.66870767) q[3];
sx q[3];
rz(0.68784586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27125204) q[0];
sx q[0];
rz(-1.7043865) q[0];
sx q[0];
rz(-1.2649076) q[0];
rz(-1.854031) q[1];
sx q[1];
rz(-0.93390211) q[1];
sx q[1];
rz(2.5087779) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.433837) q[0];
sx q[0];
rz(-2.750706) q[0];
sx q[0];
rz(0.99796064) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6980324) q[2];
sx q[2];
rz(-0.75492263) q[2];
sx q[2];
rz(-2.519543) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.22196427) q[1];
sx q[1];
rz(-1.1884513) q[1];
sx q[1];
rz(1.4098806) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0086083) q[3];
sx q[3];
rz(-0.62963943) q[3];
sx q[3];
rz(1.255577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5536993) q[2];
sx q[2];
rz(-0.18612315) q[2];
sx q[2];
rz(0.57559377) q[2];
rz(-2.0268188) q[3];
sx q[3];
rz(-1.2131178) q[3];
sx q[3];
rz(-2.7093844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1681528) q[0];
sx q[0];
rz(-1.3896717) q[0];
sx q[0];
rz(2.6218276) q[0];
rz(-1.686056) q[1];
sx q[1];
rz(-2.3958903) q[1];
sx q[1];
rz(-1.0882475) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3322743) q[0];
sx q[0];
rz(-1.921321) q[0];
sx q[0];
rz(-2.4615898) q[0];
rz(-1.9495522) q[2];
sx q[2];
rz(-2.566545) q[2];
sx q[2];
rz(1.4894007) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.714915) q[1];
sx q[1];
rz(-1.2132056) q[1];
sx q[1];
rz(-2.5044051) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6048221) q[3];
sx q[3];
rz(-1.7884975) q[3];
sx q[3];
rz(-0.061542573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9600642) q[2];
sx q[2];
rz(-2.0025415) q[2];
sx q[2];
rz(3.0354434) q[2];
rz(1.6541121) q[3];
sx q[3];
rz(-0.57523504) q[3];
sx q[3];
rz(-2.9851448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079025896) q[0];
sx q[0];
rz(-1.2747958) q[0];
sx q[0];
rz(-1.4803192) q[0];
rz(1.5638428) q[1];
sx q[1];
rz(-0.5636951) q[1];
sx q[1];
rz(0.021934358) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75087529) q[0];
sx q[0];
rz(-2.1712533) q[0];
sx q[0];
rz(-3.0200274) q[0];
rz(-pi) q[1];
rz(1.9490795) q[2];
sx q[2];
rz(-1.6334051) q[2];
sx q[2];
rz(0.68990842) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1749461) q[1];
sx q[1];
rz(-1.7632329) q[1];
sx q[1];
rz(-2.220146) q[1];
x q[2];
rz(-1.1975953) q[3];
sx q[3];
rz(-1.6587199) q[3];
sx q[3];
rz(2.4628291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.03453001) q[2];
sx q[2];
rz(-1.7835534) q[2];
sx q[2];
rz(-3.0800842) q[2];
rz(2.6574078) q[3];
sx q[3];
rz(-1.8248841) q[3];
sx q[3];
rz(2.9137602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7059962) q[0];
sx q[0];
rz(-2.0441002) q[0];
sx q[0];
rz(0.47501269) q[0];
rz(-1.8571521) q[1];
sx q[1];
rz(-0.78452763) q[1];
sx q[1];
rz(2.6433943) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9371004) q[0];
sx q[0];
rz(-1.0255359) q[0];
sx q[0];
rz(1.7496566) q[0];
rz(-pi) q[1];
x q[1];
rz(0.022749697) q[2];
sx q[2];
rz(-0.91500797) q[2];
sx q[2];
rz(-1.6237669) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.63677728) q[1];
sx q[1];
rz(-2.9801629) q[1];
sx q[1];
rz(2.982711) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3257019) q[3];
sx q[3];
rz(-2.0384841) q[3];
sx q[3];
rz(2.4809458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4755134) q[2];
sx q[2];
rz(-2.3011415) q[2];
sx q[2];
rz(-1.6125512) q[2];
rz(-2.9764002) q[3];
sx q[3];
rz(-1.8742722) q[3];
sx q[3];
rz(0.40273777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85544473) q[0];
sx q[0];
rz(-0.5583455) q[0];
sx q[0];
rz(2.1530491) q[0];
rz(-2.5033011) q[1];
sx q[1];
rz(-2.0604362) q[1];
sx q[1];
rz(0.036103006) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38449088) q[0];
sx q[0];
rz(-0.57390139) q[0];
sx q[0];
rz(1.4623653) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6240221) q[2];
sx q[2];
rz(-2.0985275) q[2];
sx q[2];
rz(-0.55014683) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7178967) q[1];
sx q[1];
rz(-1.0362018) q[1];
sx q[1];
rz(0.56294341) q[1];
rz(1.7965616) q[3];
sx q[3];
rz(-1.3610971) q[3];
sx q[3];
rz(2.5906627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4315167) q[2];
sx q[2];
rz(-2.4301131) q[2];
sx q[2];
rz(-0.21314387) q[2];
rz(1.824481) q[3];
sx q[3];
rz(-2.9340543) q[3];
sx q[3];
rz(0.801314) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31981907) q[0];
sx q[0];
rz(-1.736301) q[0];
sx q[0];
rz(2.2176493) q[0];
rz(0.79628235) q[1];
sx q[1];
rz(-2.2932107) q[1];
sx q[1];
rz(-0.36416818) q[1];
rz(-0.91847112) q[2];
sx q[2];
rz(-0.9705636) q[2];
sx q[2];
rz(2.8080429) q[2];
rz(0.65609559) q[3];
sx q[3];
rz(-0.96047209) q[3];
sx q[3];
rz(-0.44580662) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
