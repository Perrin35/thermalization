OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.035064) q[0];
sx q[0];
rz(-2.0523235) q[0];
sx q[0];
rz(2.9805592) q[0];
rz(1.610202) q[1];
sx q[1];
rz(-0.47709563) q[1];
sx q[1];
rz(0.49638003) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3408605) q[0];
sx q[0];
rz(-0.98156089) q[0];
sx q[0];
rz(0.13829921) q[0];
rz(-0.33978396) q[2];
sx q[2];
rz(-1.7614363) q[2];
sx q[2];
rz(-1.2059739) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1094184) q[1];
sx q[1];
rz(-0.86583455) q[1];
sx q[1];
rz(-1.6900307) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22827893) q[3];
sx q[3];
rz(-0.41729673) q[3];
sx q[3];
rz(-1.8018064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51497841) q[2];
sx q[2];
rz(-0.86003059) q[2];
sx q[2];
rz(-2.853945) q[2];
rz(1.7488165) q[3];
sx q[3];
rz(-2.1578622) q[3];
sx q[3];
rz(2.0387409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2709133) q[0];
sx q[0];
rz(-0.7551071) q[0];
sx q[0];
rz(2.5640008) q[0];
rz(1.5354935) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(1.5637406) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9614812) q[0];
sx q[0];
rz(-0.24678386) q[0];
sx q[0];
rz(-0.47829511) q[0];
rz(1.8446484) q[2];
sx q[2];
rz(-1.2784064) q[2];
sx q[2];
rz(-1.3016303) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2240552) q[1];
sx q[1];
rz(-1.1174669) q[1];
sx q[1];
rz(0.097150306) q[1];
x q[2];
rz(-1.413826) q[3];
sx q[3];
rz(-1.4179215) q[3];
sx q[3];
rz(2.4863941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0343798) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(2.0593026) q[2];
rz(-1.9783431) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(-2.5015586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0619693) q[0];
sx q[0];
rz(-0.29456961) q[0];
sx q[0];
rz(0.18297718) q[0];
rz(-0.043116365) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(-1.144369) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9336752) q[0];
sx q[0];
rz(-1.1753923) q[0];
sx q[0];
rz(1.0008706) q[0];
rz(-pi) q[1];
rz(1.794533) q[2];
sx q[2];
rz(-1.1615331) q[2];
sx q[2];
rz(-2.4350016) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5444762) q[1];
sx q[1];
rz(-2.5529478) q[1];
sx q[1];
rz(-2.4878923) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8539682) q[3];
sx q[3];
rz(-1.5737185) q[3];
sx q[3];
rz(0.50141108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1220876) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(0.90467492) q[2];
rz(-2.4217862) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(2.3864746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4633789) q[0];
sx q[0];
rz(-0.97020522) q[0];
sx q[0];
rz(2.4257207) q[0];
rz(0.32575682) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(-0.99266565) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2199812) q[0];
sx q[0];
rz(-2.5489759) q[0];
sx q[0];
rz(2.6556334) q[0];
rz(-0.9348346) q[2];
sx q[2];
rz(-1.0996498) q[2];
sx q[2];
rz(-2.3061894) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0996453) q[1];
sx q[1];
rz(-1.4775659) q[1];
sx q[1];
rz(-0.16361841) q[1];
rz(-pi) q[2];
rz(0.4977787) q[3];
sx q[3];
rz(-1.0504424) q[3];
sx q[3];
rz(-1.2342412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0585534) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(-1.4286208) q[2];
rz(0.84609091) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(-1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59584004) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(-2.3748421) q[0];
rz(1.2402361) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(2.7899182) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2857916) q[0];
sx q[0];
rz(-1.8579146) q[0];
sx q[0];
rz(-0.65335269) q[0];
rz(-pi) q[1];
rz(2.211116) q[2];
sx q[2];
rz(-0.74379197) q[2];
sx q[2];
rz(3.004068) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5366718) q[1];
sx q[1];
rz(-2.7465638) q[1];
sx q[1];
rz(-3.1035963) q[1];
x q[2];
rz(-1.0655754) q[3];
sx q[3];
rz(-1.8540283) q[3];
sx q[3];
rz(-0.035988228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9436283) q[2];
sx q[2];
rz(-1.6683656) q[2];
sx q[2];
rz(-2.6082805) q[2];
rz(-2.7126281) q[3];
sx q[3];
rz(-2.6140489) q[3];
sx q[3];
rz(1.126948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0241942) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(0.54164106) q[0];
rz(-2.533124) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(1.0995964) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27605155) q[0];
sx q[0];
rz(-1.9380894) q[0];
sx q[0];
rz(2.0588576) q[0];
x q[1];
rz(0.18892388) q[2];
sx q[2];
rz(-2.0940229) q[2];
sx q[2];
rz(0.97666937) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1814327) q[1];
sx q[1];
rz(-2.2391041) q[1];
sx q[1];
rz(1.798435) q[1];
x q[2];
rz(-2.3533456) q[3];
sx q[3];
rz(-1.3282093) q[3];
sx q[3];
rz(-1.6227674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6301443) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(-0.28298322) q[2];
rz(2.4354368) q[3];
sx q[3];
rz(-1.5420087) q[3];
sx q[3];
rz(2.506536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054984897) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(1.6116066) q[0];
rz(0.64487547) q[1];
sx q[1];
rz(-2.6480643) q[1];
sx q[1];
rz(2.9842916) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10282117) q[0];
sx q[0];
rz(-2.4024995) q[0];
sx q[0];
rz(-0.13939136) q[0];
x q[1];
rz(0.83432014) q[2];
sx q[2];
rz(-1.7322707) q[2];
sx q[2];
rz(1.8404567) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8017756) q[1];
sx q[1];
rz(-2.1784557) q[1];
sx q[1];
rz(-2.9031309) q[1];
x q[2];
rz(-2.5141684) q[3];
sx q[3];
rz(-1.3324696) q[3];
sx q[3];
rz(-1.081092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9592231) q[2];
sx q[2];
rz(-0.59252858) q[2];
sx q[2];
rz(-0.83089337) q[2];
rz(-0.043878555) q[3];
sx q[3];
rz(-1.9079804) q[3];
sx q[3];
rz(0.20251814) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6458994) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(-0.010852531) q[0];
rz(-0.27663484) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(-0.29327926) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54942185) q[0];
sx q[0];
rz(-2.5991419) q[0];
sx q[0];
rz(-3.1300934) q[0];
rz(-pi) q[1];
rz(-1.3515527) q[2];
sx q[2];
rz(-0.65537894) q[2];
sx q[2];
rz(0.13742451) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.20626942) q[1];
sx q[1];
rz(-0.85201293) q[1];
sx q[1];
rz(-2.1410336) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46383143) q[3];
sx q[3];
rz(-1.2004735) q[3];
sx q[3];
rz(-0.90255373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.098112) q[2];
sx q[2];
rz(-0.094823368) q[2];
sx q[2];
rz(-0.088767178) q[2];
rz(0.25012112) q[3];
sx q[3];
rz(-1.1238031) q[3];
sx q[3];
rz(-0.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1373238) q[0];
sx q[0];
rz(-0.22432888) q[0];
sx q[0];
rz(2.0776757) q[0];
rz(2.6990199) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(-1.3508505) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3063072) q[0];
sx q[0];
rz(-1.4915691) q[0];
sx q[0];
rz(-1.6552734) q[0];
rz(-pi) q[1];
rz(0.66019085) q[2];
sx q[2];
rz(-1.5056416) q[2];
sx q[2];
rz(2.8934663) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.46637725) q[1];
sx q[1];
rz(-1.5870759) q[1];
sx q[1];
rz(1.1364163) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7926932) q[3];
sx q[3];
rz(-1.4699234) q[3];
sx q[3];
rz(1.6290806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.110048) q[2];
sx q[2];
rz(-1.2908547) q[2];
sx q[2];
rz(2.6943977) q[2];
rz(1.3859008) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(2.6269004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31496012) q[0];
sx q[0];
rz(-1.5216014) q[0];
sx q[0];
rz(0.78053027) q[0];
rz(2.7138117) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(-2.9719877) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0376301) q[0];
sx q[0];
rz(-0.61015218) q[0];
sx q[0];
rz(0.9141586) q[0];
x q[1];
rz(-1.5461966) q[2];
sx q[2];
rz(-1.9925756) q[2];
sx q[2];
rz(-1.8116902) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.032220275) q[1];
sx q[1];
rz(-1.5240372) q[1];
sx q[1];
rz(2.1857775) q[1];
rz(-pi) q[2];
rz(-2.9986936) q[3];
sx q[3];
rz(-1.2522962) q[3];
sx q[3];
rz(-2.3058476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.85049373) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(0.69236857) q[2];
rz(2.678357) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(2.0326116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8511843) q[0];
sx q[0];
rz(-1.6155227) q[0];
sx q[0];
rz(0.68646705) q[0];
rz(0.51207536) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(0.26906536) q[2];
sx q[2];
rz(-1.8620327) q[2];
sx q[2];
rz(-3.0897683) q[2];
rz(-2.0797707) q[3];
sx q[3];
rz(-0.23734262) q[3];
sx q[3];
rz(1.6457641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];