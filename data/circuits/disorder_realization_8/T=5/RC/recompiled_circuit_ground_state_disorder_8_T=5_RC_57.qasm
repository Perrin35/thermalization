OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.79787624) q[0];
sx q[0];
rz(-0.98348445) q[0];
sx q[0];
rz(-0.19177344) q[0];
rz(-2.6311488) q[1];
sx q[1];
rz(-1.3671083) q[1];
sx q[1];
rz(1.3956611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53419848) q[0];
sx q[0];
rz(-1.5572131) q[0];
sx q[0];
rz(3.109038) q[0];
rz(-2.9119266) q[2];
sx q[2];
rz(-1.3896835) q[2];
sx q[2];
rz(1.0223946) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3049136) q[1];
sx q[1];
rz(-0.80614754) q[1];
sx q[1];
rz(-2.8668001) q[1];
x q[2];
rz(1.9590501) q[3];
sx q[3];
rz(-0.40273977) q[3];
sx q[3];
rz(-0.54408636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.50769606) q[2];
sx q[2];
rz(-0.32749978) q[2];
sx q[2];
rz(2.2155649) q[2];
rz(-1.5121459) q[3];
sx q[3];
rz(-2.47086) q[3];
sx q[3];
rz(-1.1431471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.73285) q[0];
sx q[0];
rz(-0.73246211) q[0];
sx q[0];
rz(0.23250411) q[0];
rz(2.0022424) q[1];
sx q[1];
rz(-1.5683697) q[1];
sx q[1];
rz(2.2154636) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3207729) q[0];
sx q[0];
rz(-2.8348587) q[0];
sx q[0];
rz(0.83267468) q[0];
rz(-0.32759362) q[2];
sx q[2];
rz(-1.6632051) q[2];
sx q[2];
rz(-1.365923) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3327712) q[1];
sx q[1];
rz(-1.0406245) q[1];
sx q[1];
rz(-1.8103241) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6517209) q[3];
sx q[3];
rz(-2.4656396) q[3];
sx q[3];
rz(2.5328758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6050379) q[2];
sx q[2];
rz(-1.522541) q[2];
sx q[2];
rz(0.79263318) q[2];
rz(0.41147453) q[3];
sx q[3];
rz(-2.1992407) q[3];
sx q[3];
rz(-2.3365432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8234392) q[0];
sx q[0];
rz(-1.0018438) q[0];
sx q[0];
rz(-0.26082984) q[0];
rz(2.8584495) q[1];
sx q[1];
rz(-1.6915551) q[1];
sx q[1];
rz(2.0859437) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2973141) q[0];
sx q[0];
rz(-1.7417272) q[0];
sx q[0];
rz(2.2173082) q[0];
x q[1];
rz(1.5745542) q[2];
sx q[2];
rz(-0.93766391) q[2];
sx q[2];
rz(2.1830677) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4312517) q[1];
sx q[1];
rz(-1.4652989) q[1];
sx q[1];
rz(2.4137817) q[1];
rz(-2.0664472) q[3];
sx q[3];
rz(-1.6248253) q[3];
sx q[3];
rz(-0.13755218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61567125) q[2];
sx q[2];
rz(-1.0893704) q[2];
sx q[2];
rz(1.2582568) q[2];
rz(-1.0387756) q[3];
sx q[3];
rz(-2.153331) q[3];
sx q[3];
rz(1.725215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1309758) q[0];
sx q[0];
rz(-1.5357786) q[0];
sx q[0];
rz(-0.11949874) q[0];
rz(-2.9833228) q[1];
sx q[1];
rz(-2.6893078) q[1];
sx q[1];
rz(1.4403042) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8058469) q[0];
sx q[0];
rz(-0.81253883) q[0];
sx q[0];
rz(1.0643908) q[0];
x q[1];
rz(1.9057299) q[2];
sx q[2];
rz(-1.4960297) q[2];
sx q[2];
rz(-0.78943816) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5938737) q[1];
sx q[1];
rz(-0.55127326) q[1];
sx q[1];
rz(-0.51993315) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0802832) q[3];
sx q[3];
rz(-1.882383) q[3];
sx q[3];
rz(-2.8924242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.67864546) q[2];
sx q[2];
rz(-3.0061649) q[2];
sx q[2];
rz(2.7079008) q[2];
rz(2.4667451) q[3];
sx q[3];
rz(-2.0516472) q[3];
sx q[3];
rz(1.1564144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40815142) q[0];
sx q[0];
rz(-1.0557405) q[0];
sx q[0];
rz(-0.59930402) q[0];
rz(-0.76404461) q[1];
sx q[1];
rz(-1.0300449) q[1];
sx q[1];
rz(-2.5291671) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9117994) q[0];
sx q[0];
rz(-0.62055885) q[0];
sx q[0];
rz(-1.3156652) q[0];
rz(-pi) q[1];
rz(0.00084288518) q[2];
sx q[2];
rz(-0.71686059) q[2];
sx q[2];
rz(-1.5886024) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80390319) q[1];
sx q[1];
rz(-0.63769996) q[1];
sx q[1];
rz(2.8578812) q[1];
x q[2];
rz(-2.0908085) q[3];
sx q[3];
rz(-2.0377167) q[3];
sx q[3];
rz(-2.0869521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6969488) q[2];
sx q[2];
rz(-2.3418043) q[2];
sx q[2];
rz(-0.4701699) q[2];
rz(0.36559513) q[3];
sx q[3];
rz(-1.9496893) q[3];
sx q[3];
rz(-2.5308334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81569165) q[0];
sx q[0];
rz(-2.3453562) q[0];
sx q[0];
rz(0.66194397) q[0];
rz(-3.0603307) q[1];
sx q[1];
rz(-2.6632402) q[1];
sx q[1];
rz(-3.001396) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7707342) q[0];
sx q[0];
rz(-2.6261289) q[0];
sx q[0];
rz(0.45592842) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1344994) q[2];
sx q[2];
rz(-2.4738174) q[2];
sx q[2];
rz(0.89950022) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0602136) q[1];
sx q[1];
rz(-1.5790931) q[1];
sx q[1];
rz(1.0992728) q[1];
rz(-2.7922157) q[3];
sx q[3];
rz(-1.3991068) q[3];
sx q[3];
rz(-2.6956345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7774272) q[2];
sx q[2];
rz(-1.3157996) q[2];
sx q[2];
rz(2.8167456) q[2];
rz(-2.9835564) q[3];
sx q[3];
rz(-0.41301781) q[3];
sx q[3];
rz(-2.1260927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82191104) q[0];
sx q[0];
rz(-3.0653636) q[0];
sx q[0];
rz(0.88777375) q[0];
rz(2.7582788) q[1];
sx q[1];
rz(-0.8822459) q[1];
sx q[1];
rz(0.15957889) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051285714) q[0];
sx q[0];
rz(-1.2538099) q[0];
sx q[0];
rz(-0.27329926) q[0];
rz(-0.55241199) q[2];
sx q[2];
rz(-1.0670976) q[2];
sx q[2];
rz(2.9744862) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7516874) q[1];
sx q[1];
rz(-1.4426822) q[1];
sx q[1];
rz(-0.37813487) q[1];
x q[2];
rz(-1.4042312) q[3];
sx q[3];
rz(-1.8565531) q[3];
sx q[3];
rz(2.8050202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54922709) q[2];
sx q[2];
rz(-1.1575907) q[2];
sx q[2];
rz(1.4166191) q[2];
rz(1.8727411) q[3];
sx q[3];
rz(-2.3609991) q[3];
sx q[3];
rz(0.95808539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45605993) q[0];
sx q[0];
rz(-2.0262418) q[0];
sx q[0];
rz(2.2400895) q[0];
rz(-2.447336) q[1];
sx q[1];
rz(-1.4638823) q[1];
sx q[1];
rz(-2.0051956) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2456348) q[0];
sx q[0];
rz(-0.16366009) q[0];
sx q[0];
rz(1.3748522) q[0];
rz(0.56635277) q[2];
sx q[2];
rz(-2.512062) q[2];
sx q[2];
rz(-0.16816631) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9688837) q[1];
sx q[1];
rz(-2.4447021) q[1];
sx q[1];
rz(-1.642176) q[1];
rz(-pi) q[2];
rz(2.091892) q[3];
sx q[3];
rz(-0.9113833) q[3];
sx q[3];
rz(-1.570861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.14094341) q[2];
sx q[2];
rz(-1.2691754) q[2];
sx q[2];
rz(-0.81306523) q[2];
rz(0.55417577) q[3];
sx q[3];
rz(-1.412609) q[3];
sx q[3];
rz(2.7798142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56162214) q[0];
sx q[0];
rz(-1.1142092) q[0];
sx q[0];
rz(2.9680874) q[0];
rz(0.42501998) q[1];
sx q[1];
rz(-1.605502) q[1];
sx q[1];
rz(0.82957155) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4886252) q[0];
sx q[0];
rz(-2.2982631) q[0];
sx q[0];
rz(-1.3759111) q[0];
x q[1];
rz(-0.046421703) q[2];
sx q[2];
rz(-1.6116465) q[2];
sx q[2];
rz(0.53120261) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8977938) q[1];
sx q[1];
rz(-1.3659371) q[1];
sx q[1];
rz(0.61989354) q[1];
x q[2];
rz(-2.9584269) q[3];
sx q[3];
rz(-0.34160994) q[3];
sx q[3];
rz(-2.6721725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7912264) q[2];
sx q[2];
rz(-1.3000877) q[2];
sx q[2];
rz(-1.680797) q[2];
rz(-2.1469877) q[3];
sx q[3];
rz(-2.1456783) q[3];
sx q[3];
rz(-2.2554956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17722546) q[0];
sx q[0];
rz(-0.57562861) q[0];
sx q[0];
rz(0.10203578) q[0];
rz(2.521934) q[1];
sx q[1];
rz(-1.6435868) q[1];
sx q[1];
rz(-0.59725753) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6031979) q[0];
sx q[0];
rz(-1.704633) q[0];
sx q[0];
rz(-2.5413245) q[0];
rz(-pi) q[1];
rz(-0.37924699) q[2];
sx q[2];
rz(-0.37988099) q[2];
sx q[2];
rz(2.3178562) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0400552) q[1];
sx q[1];
rz(-2.1360847) q[1];
sx q[1];
rz(-0.043626309) q[1];
rz(-pi) q[2];
rz(-1.6035613) q[3];
sx q[3];
rz(-2.0420342) q[3];
sx q[3];
rz(-0.12880381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20798802) q[2];
sx q[2];
rz(-2.7808069) q[2];
sx q[2];
rz(0.92850816) q[2];
rz(-1.9814631) q[3];
sx q[3];
rz(-1.7103633) q[3];
sx q[3];
rz(-0.76735705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74850294) q[0];
sx q[0];
rz(-0.7001644) q[0];
sx q[0];
rz(-2.9099332) q[0];
rz(2.0139991) q[1];
sx q[1];
rz(-0.53378202) q[1];
sx q[1];
rz(3.1376874) q[1];
rz(-2.247159) q[2];
sx q[2];
rz(-1.5147687) q[2];
sx q[2];
rz(1.1942836) q[2];
rz(-2.7990758) q[3];
sx q[3];
rz(-2.2377662) q[3];
sx q[3];
rz(-0.025442414) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
