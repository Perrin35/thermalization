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
rz(0.32349411) q[0];
sx q[0];
rz(-0.20322023) q[0];
sx q[0];
rz(0.30015552) q[0];
rz(2.7872941) q[1];
sx q[1];
rz(-1.2343255) q[1];
sx q[1];
rz(-0.77114463) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3187689) q[0];
sx q[0];
rz(-1.7576801) q[0];
sx q[0];
rz(2.9479638) q[0];
rz(-pi) q[1];
rz(-2.6080898) q[2];
sx q[2];
rz(-1.9868104) q[2];
sx q[2];
rz(-2.5737284) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1465774) q[1];
sx q[1];
rz(-1.9781917) q[1];
sx q[1];
rz(-0.40991524) q[1];
x q[2];
rz(-1.425779) q[3];
sx q[3];
rz(-2.7035664) q[3];
sx q[3];
rz(1.5256385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.99042088) q[2];
sx q[2];
rz(-2.6800818) q[2];
sx q[2];
rz(-1.0837466) q[2];
rz(-0.57461965) q[3];
sx q[3];
rz(-1.5950404) q[3];
sx q[3];
rz(-2.3641724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3144749) q[0];
sx q[0];
rz(-0.62905335) q[0];
sx q[0];
rz(2.7431059) q[0];
rz(1.9972948) q[1];
sx q[1];
rz(-0.60595787) q[1];
sx q[1];
rz(0.95432895) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0701369) q[0];
sx q[0];
rz(-3.1370039) q[0];
sx q[0];
rz(3.0194886) q[0];
rz(-pi) q[1];
rz(-1.6923087) q[2];
sx q[2];
rz(-2.4803616) q[2];
sx q[2];
rz(-1.1053156) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1283537) q[1];
sx q[1];
rz(-1.3283821) q[1];
sx q[1];
rz(-2.1552312) q[1];
rz(-0.70084533) q[3];
sx q[3];
rz(-2.4748445) q[3];
sx q[3];
rz(-1.6934044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9053456) q[2];
sx q[2];
rz(-0.50062847) q[2];
sx q[2];
rz(-0.11032571) q[2];
rz(-2.3162383) q[3];
sx q[3];
rz(-0.76485157) q[3];
sx q[3];
rz(2.4013405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.1684619) q[0];
sx q[0];
rz(-2.923712) q[0];
sx q[0];
rz(2.2595898) q[0];
rz(-1.0285671) q[1];
sx q[1];
rz(-0.55362892) q[1];
sx q[1];
rz(-2.2291768) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19174978) q[0];
sx q[0];
rz(-2.0217359) q[0];
sx q[0];
rz(-3.0613857) q[0];
rz(-0.71850168) q[2];
sx q[2];
rz(-1.6616491) q[2];
sx q[2];
rz(-1.529734) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0076651) q[1];
sx q[1];
rz(-1.0093736) q[1];
sx q[1];
rz(-0.39427653) q[1];
rz(-pi) q[2];
rz(-2.8863593) q[3];
sx q[3];
rz(-1.1369656) q[3];
sx q[3];
rz(1.1800769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2622111) q[2];
sx q[2];
rz(-2.8935581) q[2];
sx q[2];
rz(-2.3233419) q[2];
rz(-2.7815871) q[3];
sx q[3];
rz(-1.2986978) q[3];
sx q[3];
rz(3.103783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4358732) q[0];
sx q[0];
rz(-0.95054764) q[0];
sx q[0];
rz(-1.1608423) q[0];
rz(0.97087651) q[1];
sx q[1];
rz(-2.2258046) q[1];
sx q[1];
rz(0.11347778) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9805841) q[0];
sx q[0];
rz(-0.9763813) q[0];
sx q[0];
rz(-1.0779146) q[0];
rz(-pi) q[1];
rz(1.6062097) q[2];
sx q[2];
rz(-1.5839911) q[2];
sx q[2];
rz(-0.21617568) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.66333713) q[1];
sx q[1];
rz(-2.2436444) q[1];
sx q[1];
rz(-0.70391432) q[1];
rz(0.63469074) q[3];
sx q[3];
rz(-2.1116167) q[3];
sx q[3];
rz(-1.5405003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.11161441) q[2];
sx q[2];
rz(-1.8054211) q[2];
sx q[2];
rz(0.90799904) q[2];
rz(0.32634431) q[3];
sx q[3];
rz(-0.55485266) q[3];
sx q[3];
rz(-0.42181304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.570785) q[0];
sx q[0];
rz(-2.3796005) q[0];
sx q[0];
rz(-2.1245891) q[0];
rz(-0.48655888) q[1];
sx q[1];
rz(-1.0252955) q[1];
sx q[1];
rz(0.09495458) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21790378) q[0];
sx q[0];
rz(-0.14553864) q[0];
sx q[0];
rz(0.83011629) q[0];
rz(-pi) q[1];
rz(-1.3516897) q[2];
sx q[2];
rz(-2.3833102) q[2];
sx q[2];
rz(-2.3996283) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9776756) q[1];
sx q[1];
rz(-1.9174121) q[1];
sx q[1];
rz(-1.279128) q[1];
rz(-pi) q[2];
rz(2.8311555) q[3];
sx q[3];
rz(-1.8707898) q[3];
sx q[3];
rz(-2.7315549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2056556) q[2];
sx q[2];
rz(-0.63465261) q[2];
sx q[2];
rz(2.671799) q[2];
rz(0.69822407) q[3];
sx q[3];
rz(-1.2263115) q[3];
sx q[3];
rz(-2.8214112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53223377) q[0];
sx q[0];
rz(-2.4025752) q[0];
sx q[0];
rz(0.21482378) q[0];
rz(0.23669446) q[1];
sx q[1];
rz(-2.5645945) q[1];
sx q[1];
rz(-0.41928852) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.740991) q[0];
sx q[0];
rz(-3.0508578) q[0];
sx q[0];
rz(2.0965133) q[0];
rz(-2.8307011) q[2];
sx q[2];
rz(-1.0513554) q[2];
sx q[2];
rz(-2.3497557) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9671075) q[1];
sx q[1];
rz(-1.9000783) q[1];
sx q[1];
rz(-0.21547538) q[1];
rz(-pi) q[2];
rz(1.3724519) q[3];
sx q[3];
rz(-2.2738278) q[3];
sx q[3];
rz(-1.3217712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.145891) q[2];
sx q[2];
rz(-0.13549165) q[2];
sx q[2];
rz(-3.0502012) q[2];
rz(0.35178301) q[3];
sx q[3];
rz(-2.3942949) q[3];
sx q[3];
rz(2.6763797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9212937) q[0];
sx q[0];
rz(-1.1722246) q[0];
sx q[0];
rz(-1.175513) q[0];
rz(0.57918864) q[1];
sx q[1];
rz(-2.3349031) q[1];
sx q[1];
rz(0.18956345) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3407758) q[0];
sx q[0];
rz(-0.22332668) q[0];
sx q[0];
rz(0.45036611) q[0];
x q[1];
rz(1.4681508) q[2];
sx q[2];
rz(-2.4375439) q[2];
sx q[2];
rz(-1.8542022) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.70041839) q[1];
sx q[1];
rz(-0.75276276) q[1];
sx q[1];
rz(-1.5819248) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.060154288) q[3];
sx q[3];
rz(-1.0511479) q[3];
sx q[3];
rz(2.1148256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.2438948) q[2];
sx q[2];
rz(-1.2995259) q[2];
sx q[2];
rz(-2.6356836) q[2];
rz(2.3627152) q[3];
sx q[3];
rz(-2.7454822) q[3];
sx q[3];
rz(-1.2946543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.8706354) q[0];
sx q[0];
rz(-1.0519692) q[0];
sx q[0];
rz(-1.1554931) q[0];
rz(0.0011860154) q[1];
sx q[1];
rz(-2.3540034) q[1];
sx q[1];
rz(-0.40518618) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9987273) q[0];
sx q[0];
rz(-2.0901276) q[0];
sx q[0];
rz(-1.0267535) q[0];
rz(2.9878163) q[2];
sx q[2];
rz(-1.1660053) q[2];
sx q[2];
rz(-2.5396233) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2470513) q[1];
sx q[1];
rz(-0.76746583) q[1];
sx q[1];
rz(1.6805499) q[1];
x q[2];
rz(-2.6514945) q[3];
sx q[3];
rz(-2.9967272) q[3];
sx q[3];
rz(-0.22508301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.81795168) q[2];
sx q[2];
rz(-1.1672856) q[2];
sx q[2];
rz(-0.52158588) q[2];
rz(-3.0449872) q[3];
sx q[3];
rz(-2.1229027) q[3];
sx q[3];
rz(-0.49083138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38799649) q[0];
sx q[0];
rz(-2.2612408) q[0];
sx q[0];
rz(-3.0269347) q[0];
rz(-1.8278587) q[1];
sx q[1];
rz(-1.4559682) q[1];
sx q[1];
rz(-0.1350666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5324482) q[0];
sx q[0];
rz(-1.8290863) q[0];
sx q[0];
rz(-1.7946971) q[0];
x q[1];
rz(1.8377322) q[2];
sx q[2];
rz(-1.9501172) q[2];
sx q[2];
rz(2.1960559) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5852063) q[1];
sx q[1];
rz(-1.408768) q[1];
sx q[1];
rz(-0.94604793) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6879795) q[3];
sx q[3];
rz(-2.5796267) q[3];
sx q[3];
rz(-2.9148852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.31585285) q[2];
sx q[2];
rz(-1.9261381) q[2];
sx q[2];
rz(2.1960171) q[2];
rz(-0.25997508) q[3];
sx q[3];
rz(-2.1187449) q[3];
sx q[3];
rz(-1.1280577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.05397) q[0];
sx q[0];
rz(-1.0902417) q[0];
sx q[0];
rz(-1.7386275) q[0];
rz(-0.91579092) q[1];
sx q[1];
rz(-2.5251838) q[1];
sx q[1];
rz(-0.91555196) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4175891) q[0];
sx q[0];
rz(-1.4737098) q[0];
sx q[0];
rz(-1.7018945) q[0];
rz(2.0924641) q[2];
sx q[2];
rz(-1.7384439) q[2];
sx q[2];
rz(1.1225134) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7448235) q[1];
sx q[1];
rz(-1.1698548) q[1];
sx q[1];
rz(-2.045782) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9718982) q[3];
sx q[3];
rz(-0.86830074) q[3];
sx q[3];
rz(-2.9164721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0182858) q[2];
sx q[2];
rz(-2.4985963) q[2];
sx q[2];
rz(-0.40261889) q[2];
rz(1.1766524) q[3];
sx q[3];
rz(-0.28665701) q[3];
sx q[3];
rz(-2.3736242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1422414) q[0];
sx q[0];
rz(-0.9724697) q[0];
sx q[0];
rz(-1.019626) q[0];
rz(-0.6877407) q[1];
sx q[1];
rz(-2.3383457) q[1];
sx q[1];
rz(1.8170423) q[1];
rz(-2.4551212) q[2];
sx q[2];
rz(-2.1778637) q[2];
sx q[2];
rz(-1.7354497) q[2];
rz(-2.5108093) q[3];
sx q[3];
rz(-1.3628241) q[3];
sx q[3];
rz(0.48490111) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
