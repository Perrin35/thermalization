OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82792038) q[0];
sx q[0];
rz(1.496324) q[0];
sx q[0];
rz(9.6376251) q[0];
rz(1.7362453) q[1];
sx q[1];
rz(4.6680556) q[1];
sx q[1];
rz(13.144346) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62993193) q[0];
sx q[0];
rz(-1.1714456) q[0];
sx q[0];
rz(2.6005247) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18696286) q[2];
sx q[2];
rz(-2.2622795) q[2];
sx q[2];
rz(2.5656467) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.69691519) q[1];
sx q[1];
rz(-1.5648876) q[1];
sx q[1];
rz(1.1586994) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37761117) q[3];
sx q[3];
rz(-2.7481573) q[3];
sx q[3];
rz(-2.4675219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9389682) q[2];
sx q[2];
rz(-1.6309996) q[2];
sx q[2];
rz(-2.5228339) q[2];
rz(-2.7135811) q[3];
sx q[3];
rz(-1.908554) q[3];
sx q[3];
rz(1.4324043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.0851058) q[0];
sx q[0];
rz(-0.036186941) q[0];
sx q[0];
rz(2.4888743) q[0];
rz(-2.808049) q[1];
sx q[1];
rz(-2.3338552) q[1];
sx q[1];
rz(-0.52887598) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31007209) q[0];
sx q[0];
rz(-1.0640826) q[0];
sx q[0];
rz(-2.7181198) q[0];
rz(1.7414167) q[2];
sx q[2];
rz(-0.58310027) q[2];
sx q[2];
rz(1.0379593) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6436504) q[1];
sx q[1];
rz(-1.4406246) q[1];
sx q[1];
rz(1.4243717) q[1];
rz(-0.06270285) q[3];
sx q[3];
rz(-2.3093105) q[3];
sx q[3];
rz(0.073918561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7876106) q[2];
sx q[2];
rz(-1.748961) q[2];
sx q[2];
rz(-0.37484136) q[2];
rz(-1.9311054) q[3];
sx q[3];
rz(-2.6755302) q[3];
sx q[3];
rz(1.9207641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1061363) q[0];
sx q[0];
rz(-1.9566256) q[0];
sx q[0];
rz(0.55738604) q[0];
rz(-2.4222971) q[1];
sx q[1];
rz(-2.6917515) q[1];
sx q[1];
rz(0.61190277) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.681088) q[0];
sx q[0];
rz(-1.2375481) q[0];
sx q[0];
rz(-1.2171576) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8863897) q[2];
sx q[2];
rz(-1.4006613) q[2];
sx q[2];
rz(-1.1332133) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0825838) q[1];
sx q[1];
rz(-1.6477668) q[1];
sx q[1];
rz(0.99748912) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7313084) q[3];
sx q[3];
rz(-0.43528523) q[3];
sx q[3];
rz(1.9225863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3376075) q[2];
sx q[2];
rz(-1.7943725) q[2];
sx q[2];
rz(-0.95532974) q[2];
rz(0.56940007) q[3];
sx q[3];
rz(-2.0208713) q[3];
sx q[3];
rz(1.7795405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78205458) q[0];
sx q[0];
rz(-2.8466917) q[0];
sx q[0];
rz(-3.0468347) q[0];
rz(-1.5358198) q[1];
sx q[1];
rz(-1.1631807) q[1];
sx q[1];
rz(-2.7045429) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8931797) q[0];
sx q[0];
rz(-1.4709269) q[0];
sx q[0];
rz(3.0573387) q[0];
rz(-pi) q[1];
rz(-2.7669698) q[2];
sx q[2];
rz(-2.6772876) q[2];
sx q[2];
rz(1.2863867) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.028965) q[1];
sx q[1];
rz(-1.3092759) q[1];
sx q[1];
rz(-3.1413743) q[1];
x q[2];
rz(1.0867723) q[3];
sx q[3];
rz(-0.71992249) q[3];
sx q[3];
rz(0.37527028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7729418) q[2];
sx q[2];
rz(-2.3029885) q[2];
sx q[2];
rz(1.8013901) q[2];
rz(-2.7784427) q[3];
sx q[3];
rz(-0.64195091) q[3];
sx q[3];
rz(0.44410458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3602314) q[0];
sx q[0];
rz(-2.4960127) q[0];
sx q[0];
rz(0.050882291) q[0];
rz(2.6119192) q[1];
sx q[1];
rz(-2.5310204) q[1];
sx q[1];
rz(-3.1231336) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0063361703) q[0];
sx q[0];
rz(-2.0783892) q[0];
sx q[0];
rz(-2.7897842) q[0];
x q[1];
rz(-1.7713806) q[2];
sx q[2];
rz(-2.2772352) q[2];
sx q[2];
rz(0.7274703) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6405271) q[1];
sx q[1];
rz(-2.200715) q[1];
sx q[1];
rz(-0.34754158) q[1];
rz(-pi) q[2];
rz(2.6191447) q[3];
sx q[3];
rz(-1.7237437) q[3];
sx q[3];
rz(-2.3207823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9336443) q[2];
sx q[2];
rz(-1.4241781) q[2];
sx q[2];
rz(-1.8682182) q[2];
rz(-0.4452855) q[3];
sx q[3];
rz(-0.71355009) q[3];
sx q[3];
rz(-0.85550296) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91525602) q[0];
sx q[0];
rz(-0.38050088) q[0];
sx q[0];
rz(-2.8476727) q[0];
rz(-1.8982915) q[1];
sx q[1];
rz(-1.2970122) q[1];
sx q[1];
rz(-1.6506724) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3476171) q[0];
sx q[0];
rz(-1.7375943) q[0];
sx q[0];
rz(2.1991437) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4751728) q[2];
sx q[2];
rz(-1.2461414) q[2];
sx q[2];
rz(2.4839475) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8141802) q[1];
sx q[1];
rz(-2.8490863) q[1];
sx q[1];
rz(2.4027225) q[1];
rz(-pi) q[2];
rz(0.22533234) q[3];
sx q[3];
rz(-1.1374416) q[3];
sx q[3];
rz(-3.0399377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2095498) q[2];
sx q[2];
rz(-1.2806634) q[2];
sx q[2];
rz(0.10188421) q[2];
rz(2.3515676) q[3];
sx q[3];
rz(-0.70370379) q[3];
sx q[3];
rz(1.8858645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1387966) q[0];
sx q[0];
rz(-0.094745435) q[0];
sx q[0];
rz(-2.5139659) q[0];
rz(1.3603285) q[1];
sx q[1];
rz(-1.5831213) q[1];
sx q[1];
rz(0.21242177) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4110574) q[0];
sx q[0];
rz(-0.90170398) q[0];
sx q[0];
rz(-0.47242609) q[0];
rz(0.97456206) q[2];
sx q[2];
rz(-1.803054) q[2];
sx q[2];
rz(-1.5982826) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0052988) q[1];
sx q[1];
rz(-1.7583656) q[1];
sx q[1];
rz(-1.0636399) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3644629) q[3];
sx q[3];
rz(-0.6327968) q[3];
sx q[3];
rz(1.5416886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3820496) q[2];
sx q[2];
rz(-2.2727649) q[2];
sx q[2];
rz(0.36054912) q[2];
rz(0.77066317) q[3];
sx q[3];
rz(-1.2040851) q[3];
sx q[3];
rz(0.26309052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78779546) q[0];
sx q[0];
rz(-2.0833092) q[0];
sx q[0];
rz(-1.298792) q[0];
rz(-2.6020488) q[1];
sx q[1];
rz(-1.6863457) q[1];
sx q[1];
rz(-1.4849327) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.36168) q[0];
sx q[0];
rz(-0.64881262) q[0];
sx q[0];
rz(-0.54260962) q[0];
rz(2.7168112) q[2];
sx q[2];
rz(-1.3662587) q[2];
sx q[2];
rz(-1.0077602) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6816652) q[1];
sx q[1];
rz(-1.2673389) q[1];
sx q[1];
rz(-3.1292679) q[1];
rz(-pi) q[2];
rz(0.45739037) q[3];
sx q[3];
rz(-0.76888385) q[3];
sx q[3];
rz(-0.035280003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7654884) q[2];
sx q[2];
rz(-1.2771983) q[2];
sx q[2];
rz(-0.45846024) q[2];
rz(-1.9930528) q[3];
sx q[3];
rz(-0.12980041) q[3];
sx q[3];
rz(0.54500088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.5371567) q[0];
sx q[0];
rz(-0.24188365) q[0];
sx q[0];
rz(-2.3271374) q[0];
rz(0.7704598) q[1];
sx q[1];
rz(-1.6043681) q[1];
sx q[1];
rz(1.8155712) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4453989) q[0];
sx q[0];
rz(-2.8210381) q[0];
sx q[0];
rz(0.93648367) q[0];
rz(-3.0203462) q[2];
sx q[2];
rz(-0.7077982) q[2];
sx q[2];
rz(1.540465) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0596366) q[1];
sx q[1];
rz(-0.15742144) q[1];
sx q[1];
rz(-1.5035968) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3877478) q[3];
sx q[3];
rz(-1.7244851) q[3];
sx q[3];
rz(2.1240687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.068437964) q[2];
sx q[2];
rz(-0.96403733) q[2];
sx q[2];
rz(0.55802074) q[2];
rz(-2.8716715) q[3];
sx q[3];
rz(-2.8840265) q[3];
sx q[3];
rz(0.71815467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.54616) q[0];
sx q[0];
rz(-1.8159001) q[0];
sx q[0];
rz(-2.2823855) q[0];
rz(-2.4190306) q[1];
sx q[1];
rz(-0.95046202) q[1];
sx q[1];
rz(1.7019466) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7116657) q[0];
sx q[0];
rz(-1.5646805) q[0];
sx q[0];
rz(-1.5625009) q[0];
rz(-pi) q[1];
rz(1.7990803) q[2];
sx q[2];
rz(-1.205027) q[2];
sx q[2];
rz(-1.7342784) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1240012) q[1];
sx q[1];
rz(-2.3307226) q[1];
sx q[1];
rz(1.8161185) q[1];
x q[2];
rz(-1.3297746) q[3];
sx q[3];
rz(-1.1794349) q[3];
sx q[3];
rz(0.64590871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47448164) q[2];
sx q[2];
rz(-1.2789395) q[2];
sx q[2];
rz(-1.6472335) q[2];
rz(0.4161559) q[3];
sx q[3];
rz(-1.4978724) q[3];
sx q[3];
rz(-0.070629899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011099747) q[0];
sx q[0];
rz(-2.6051705) q[0];
sx q[0];
rz(2.312533) q[0];
rz(0.88498712) q[1];
sx q[1];
rz(-2.5288455) q[1];
sx q[1];
rz(1.7987342) q[1];
rz(1.5985684) q[2];
sx q[2];
rz(-1.0226316) q[2];
sx q[2];
rz(1.2050261) q[2];
rz(1.8237413) q[3];
sx q[3];
rz(-0.80251903) q[3];
sx q[3];
rz(2.4783721) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
