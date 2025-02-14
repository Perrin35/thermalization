OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(5.7558007) q[0];
sx q[0];
rz(5.7100073) q[0];
sx q[0];
rz(10.041458) q[0];
rz(2.1332027) q[1];
sx q[1];
rz(-0.73928666) q[1];
sx q[1];
rz(0.33831212) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3544681) q[0];
sx q[0];
rz(-1.6068874) q[0];
sx q[0];
rz(-2.8990394) q[0];
rz(-2.7213547) q[2];
sx q[2];
rz(-2.3048688) q[2];
sx q[2];
rz(-0.24009304) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1627312) q[1];
sx q[1];
rz(-0.56038364) q[1];
sx q[1];
rz(2.8864229) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7565207) q[3];
sx q[3];
rz(-0.76558569) q[3];
sx q[3];
rz(-3.035745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3931291) q[2];
sx q[2];
rz(-1.7270361) q[2];
sx q[2];
rz(-0.65790042) q[2];
rz(-0.45514485) q[3];
sx q[3];
rz(-2.9908266) q[3];
sx q[3];
rz(1.8337102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2748134) q[0];
sx q[0];
rz(-1.7300737) q[0];
sx q[0];
rz(2.859512) q[0];
rz(2.8981949) q[1];
sx q[1];
rz(-1.2671821) q[1];
sx q[1];
rz(-0.74877053) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3385612) q[0];
sx q[0];
rz(-0.73006064) q[0];
sx q[0];
rz(2.8499313) q[0];
rz(-pi) q[1];
rz(1.6328598) q[2];
sx q[2];
rz(-0.12189874) q[2];
sx q[2];
rz(2.1699407) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6525014) q[1];
sx q[1];
rz(-1.6630807) q[1];
sx q[1];
rz(2.0745866) q[1];
rz(-pi) q[2];
rz(0.19939662) q[3];
sx q[3];
rz(-1.6391338) q[3];
sx q[3];
rz(-1.133356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.28027174) q[2];
sx q[2];
rz(-1.5254285) q[2];
sx q[2];
rz(1.3410428) q[2];
rz(-1.1567814) q[3];
sx q[3];
rz(-0.78514922) q[3];
sx q[3];
rz(0.7635428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84083104) q[0];
sx q[0];
rz(-0.16591993) q[0];
sx q[0];
rz(0.65735835) q[0];
rz(-1.9961458) q[1];
sx q[1];
rz(-2.1871388) q[1];
sx q[1];
rz(-2.3232536) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20387801) q[0];
sx q[0];
rz(-0.42417198) q[0];
sx q[0];
rz(0.61221497) q[0];
rz(-1.7162227) q[2];
sx q[2];
rz(-0.31842318) q[2];
sx q[2];
rz(-0.91998902) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0967674) q[1];
sx q[1];
rz(-0.92508537) q[1];
sx q[1];
rz(2.195408) q[1];
rz(-pi) q[2];
rz(-2.0243689) q[3];
sx q[3];
rz(-1.6370341) q[3];
sx q[3];
rz(0.57285786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1546617) q[2];
sx q[2];
rz(-1.8795452) q[2];
sx q[2];
rz(1.3531125) q[2];
rz(-0.96495572) q[3];
sx q[3];
rz(-1.8392287) q[3];
sx q[3];
rz(-0.42199782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44593909) q[0];
sx q[0];
rz(-2.4762479) q[0];
sx q[0];
rz(1.2009784) q[0];
rz(-3.1383842) q[1];
sx q[1];
rz(-1.708958) q[1];
sx q[1];
rz(2.9046955) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99588441) q[0];
sx q[0];
rz(-2.0170324) q[0];
sx q[0];
rz(0.97490411) q[0];
rz(0.17572524) q[2];
sx q[2];
rz(-0.64470664) q[2];
sx q[2];
rz(0.75632655) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2806816) q[1];
sx q[1];
rz(-2.1691285) q[1];
sx q[1];
rz(-2.9053238) q[1];
rz(-pi) q[2];
rz(2.3637216) q[3];
sx q[3];
rz(-2.155455) q[3];
sx q[3];
rz(-2.8608222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0995471) q[2];
sx q[2];
rz(-1.2914265) q[2];
sx q[2];
rz(0.42638865) q[2];
rz(-1.03553) q[3];
sx q[3];
rz(-1.6798881) q[3];
sx q[3];
rz(-2.5199913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7972888) q[0];
sx q[0];
rz(-3.099589) q[0];
sx q[0];
rz(2.7916743) q[0];
rz(1.2087076) q[1];
sx q[1];
rz(-1.8588926) q[1];
sx q[1];
rz(3.1399609) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37536538) q[0];
sx q[0];
rz(-1.3266801) q[0];
sx q[0];
rz(1.3993511) q[0];
rz(-pi) q[1];
rz(1.6095743) q[2];
sx q[2];
rz(-1.2735521) q[2];
sx q[2];
rz(-0.69678885) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7824911) q[1];
sx q[1];
rz(-1.3939855) q[1];
sx q[1];
rz(-3.0935118) q[1];
x q[2];
rz(-1.8408082) q[3];
sx q[3];
rz(-1.7740267) q[3];
sx q[3];
rz(-3.0940987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95167595) q[2];
sx q[2];
rz(-0.92477551) q[2];
sx q[2];
rz(0.16119371) q[2];
rz(-0.4392043) q[3];
sx q[3];
rz(-0.74857155) q[3];
sx q[3];
rz(-0.85328931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8301903) q[0];
sx q[0];
rz(-1.7364194) q[0];
sx q[0];
rz(-2.3468974) q[0];
rz(-1.7757802) q[1];
sx q[1];
rz(-2.3410485) q[1];
sx q[1];
rz(-2.6672003) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6013941) q[0];
sx q[0];
rz(-2.0776333) q[0];
sx q[0];
rz(0.43497928) q[0];
x q[1];
rz(0.55616711) q[2];
sx q[2];
rz(-0.63236299) q[2];
sx q[2];
rz(1.4390505) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5582592) q[1];
sx q[1];
rz(-0.72860241) q[1];
sx q[1];
rz(2.5069951) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9871363) q[3];
sx q[3];
rz(-1.8963834) q[3];
sx q[3];
rz(1.9501571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.61364335) q[2];
sx q[2];
rz(-1.2443292) q[2];
sx q[2];
rz(-2.6141686) q[2];
rz(-2.534965) q[3];
sx q[3];
rz(-0.18692034) q[3];
sx q[3];
rz(-2.0691779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2622751) q[0];
sx q[0];
rz(-2.9476705) q[0];
sx q[0];
rz(2.344017) q[0];
rz(-0.89353117) q[1];
sx q[1];
rz(-2.306566) q[1];
sx q[1];
rz(-2.8614047) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3994795) q[0];
sx q[0];
rz(-2.7055938) q[0];
sx q[0];
rz(1.9943555) q[0];
rz(-2.1970388) q[2];
sx q[2];
rz(-2.124386) q[2];
sx q[2];
rz(-1.9596726) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.97238648) q[1];
sx q[1];
rz(-1.5833588) q[1];
sx q[1];
rz(-1.088284) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7817621) q[3];
sx q[3];
rz(-2.1323279) q[3];
sx q[3];
rz(-0.4175182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7402652) q[2];
sx q[2];
rz(-1.836931) q[2];
sx q[2];
rz(1.2124445) q[2];
rz(0.23981833) q[3];
sx q[3];
rz(-0.66767728) q[3];
sx q[3];
rz(2.5734606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78385335) q[0];
sx q[0];
rz(-1.724406) q[0];
sx q[0];
rz(1.8200112) q[0];
rz(-0.18121885) q[1];
sx q[1];
rz(-1.8099338) q[1];
sx q[1];
rz(-0.063057335) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1940799) q[0];
sx q[0];
rz(-0.79160684) q[0];
sx q[0];
rz(-2.9400154) q[0];
rz(-pi) q[1];
rz(-1.2853773) q[2];
sx q[2];
rz(-1.6550502) q[2];
sx q[2];
rz(-0.38049305) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6136421) q[1];
sx q[1];
rz(-2.7220352) q[1];
sx q[1];
rz(-2.502524) q[1];
x q[2];
rz(-1.5884871) q[3];
sx q[3];
rz(-2.9790386) q[3];
sx q[3];
rz(-1.7944195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7276089) q[2];
sx q[2];
rz(-1.0980282) q[2];
sx q[2];
rz(-1.46924) q[2];
rz(-1.7478878) q[3];
sx q[3];
rz(-1.7017476) q[3];
sx q[3];
rz(0.20974717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53836981) q[0];
sx q[0];
rz(-2.5312238) q[0];
sx q[0];
rz(-2.1446153) q[0];
rz(-2.3485377) q[1];
sx q[1];
rz(-1.7231562) q[1];
sx q[1];
rz(-1.1096032) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1556517) q[0];
sx q[0];
rz(-1.2796667) q[0];
sx q[0];
rz(-2.554214) q[0];
rz(-2.419554) q[2];
sx q[2];
rz(-1.8991611) q[2];
sx q[2];
rz(-1.9313081) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.058179) q[1];
sx q[1];
rz(-1.021487) q[1];
sx q[1];
rz(2.5843402) q[1];
rz(-pi) q[2];
x q[2];
rz(2.505198) q[3];
sx q[3];
rz(-1.706746) q[3];
sx q[3];
rz(0.037347945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7290466) q[2];
sx q[2];
rz(-1.46526) q[2];
sx q[2];
rz(1.2726146) q[2];
rz(-1.1254958) q[3];
sx q[3];
rz(-2.9477305) q[3];
sx q[3];
rz(3.061749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7703055) q[0];
sx q[0];
rz(-1.1530387) q[0];
sx q[0];
rz(-0.09819296) q[0];
rz(-1.498361) q[1];
sx q[1];
rz(-1.1474835) q[1];
sx q[1];
rz(2.3032545) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9591374) q[0];
sx q[0];
rz(-2.3172624) q[0];
sx q[0];
rz(1.3289544) q[0];
x q[1];
rz(-0.25213253) q[2];
sx q[2];
rz(-1.6246535) q[2];
sx q[2];
rz(-2.5359652) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3915836) q[1];
sx q[1];
rz(-2.6152088) q[1];
sx q[1];
rz(2.1571804) q[1];
rz(0.43180744) q[3];
sx q[3];
rz(-1.3696967) q[3];
sx q[3];
rz(-0.74936282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7958293) q[2];
sx q[2];
rz(-0.886262) q[2];
sx q[2];
rz(0.32540992) q[2];
rz(0.77643967) q[3];
sx q[3];
rz(-2.1263945) q[3];
sx q[3];
rz(0.6071035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0008739) q[0];
sx q[0];
rz(-1.2043395) q[0];
sx q[0];
rz(2.1111063) q[0];
rz(-0.28618947) q[1];
sx q[1];
rz(-1.9356526) q[1];
sx q[1];
rz(-2.0331358) q[1];
rz(0.76079615) q[2];
sx q[2];
rz(-1.8626871) q[2];
sx q[2];
rz(-2.5260851) q[2];
rz(0.48578942) q[3];
sx q[3];
rz(-0.64668568) q[3];
sx q[3];
rz(0.034737094) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
