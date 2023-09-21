OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10652868) q[0];
sx q[0];
rz(-1.0892692) q[0];
sx q[0];
rz(-2.9805592) q[0];
rz(1.610202) q[1];
sx q[1];
rz(-0.47709563) q[1];
sx q[1];
rz(0.49638003) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8344603) q[0];
sx q[0];
rz(-1.4559329) q[0];
sx q[0];
rz(-0.97712028) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6167294) q[2];
sx q[2];
rz(-0.3877936) q[2];
sx q[2];
rz(0.85688574) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6160994) q[1];
sx q[1];
rz(-1.6615189) q[1];
sx q[1];
rz(-2.4331122) q[1];
x q[2];
rz(-1.4708038) q[3];
sx q[3];
rz(-1.1649719) q[3];
sx q[3];
rz(2.0506746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6266142) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(-2.853945) q[2];
rz(1.3927762) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(-1.1028517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2709133) q[0];
sx q[0];
rz(-0.7551071) q[0];
sx q[0];
rz(-2.5640008) q[0];
rz(-1.6060991) q[1];
sx q[1];
rz(-2.0237193) q[1];
sx q[1];
rz(-1.5637406) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9614812) q[0];
sx q[0];
rz(-2.8948088) q[0];
sx q[0];
rz(0.47829511) q[0];
x q[1];
rz(2.8385542) q[2];
sx q[2];
rz(-1.3088471) q[2];
sx q[2];
rz(0.18837243) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2240552) q[1];
sx q[1];
rz(-2.0241258) q[1];
sx q[1];
rz(3.0444423) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7277667) q[3];
sx q[3];
rz(-1.4179215) q[3];
sx q[3];
rz(0.65519858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0343798) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(1.0822901) q[2];
rz(1.9783431) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(2.5015586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0796233) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(-2.9586155) q[0];
rz(3.0984763) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(-1.144369) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2079175) q[0];
sx q[0];
rz(-1.9662004) q[0];
sx q[0];
rz(2.140722) q[0];
rz(-1.3470596) q[2];
sx q[2];
rz(-1.9800595) q[2];
sx q[2];
rz(2.4350016) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.59359264) q[1];
sx q[1];
rz(-1.9152194) q[1];
sx q[1];
rz(2.6542632) q[1];
rz(-pi) q[2];
rz(-0.28762443) q[3];
sx q[3];
rz(-1.5678741) q[3];
sx q[3];
rz(-2.6401816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1220876) q[2];
sx q[2];
rz(-1.1818037) q[2];
sx q[2];
rz(-0.90467492) q[2];
rz(-2.4217862) q[3];
sx q[3];
rz(-2.7147839) q[3];
sx q[3];
rz(-2.3864746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6782137) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(2.4257207) q[0];
rz(-0.32575682) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(-0.99266565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7870165) q[0];
sx q[0];
rz(-2.0873318) q[0];
sx q[0];
rz(-1.2660962) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2797212) q[2];
sx q[2];
rz(-2.3700691) q[2];
sx q[2];
rz(0.18415235) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0419473) q[1];
sx q[1];
rz(-1.4775659) q[1];
sx q[1];
rz(-0.16361841) q[1];
x q[2];
rz(0.87611115) q[3];
sx q[3];
rz(-2.4377341) q[3];
sx q[3];
rz(-1.0775523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0830393) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(-1.7129718) q[2];
rz(2.2955017) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(1.2231474) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5457526) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(-2.3748421) q[0];
rz(-1.2402361) q[1];
sx q[1];
rz(-0.84754544) q[1];
sx q[1];
rz(-0.3516745) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36067097) q[0];
sx q[0];
rz(-2.43649) q[0];
sx q[0];
rz(2.6893925) q[0];
x q[1];
rz(2.211116) q[2];
sx q[2];
rz(-2.3978007) q[2];
sx q[2];
rz(-3.004068) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5778351) q[1];
sx q[1];
rz(-1.9655242) q[1];
sx q[1];
rz(-1.5866336) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0655754) q[3];
sx q[3];
rz(-1.2875644) q[3];
sx q[3];
rz(3.1056044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19796431) q[2];
sx q[2];
rz(-1.6683656) q[2];
sx q[2];
rz(2.6082805) q[2];
rz(-2.7126281) q[3];
sx q[3];
rz(-2.6140489) q[3];
sx q[3];
rz(1.126948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0241942) q[0];
sx q[0];
rz(-2.0485931) q[0];
sx q[0];
rz(-2.5999516) q[0];
rz(2.533124) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(-1.0995964) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8894419) q[0];
sx q[0];
rz(-0.60177207) q[0];
sx q[0];
rz(-0.88366951) q[0];
rz(0.18892388) q[2];
sx q[2];
rz(-2.0940229) q[2];
sx q[2];
rz(-2.1649233) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.53193608) q[1];
sx q[1];
rz(-1.3927288) q[1];
sx q[1];
rz(-2.460536) q[1];
x q[2];
rz(-1.908329) q[3];
sx q[3];
rz(-2.3300155) q[3];
sx q[3];
rz(0.28901643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5114484) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(0.28298322) q[2];
rz(-2.4354368) q[3];
sx q[3];
rz(-1.5420087) q[3];
sx q[3];
rz(-2.506536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054984897) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(1.6116066) q[0];
rz(-2.4967172) q[1];
sx q[1];
rz(-2.6480643) q[1];
sx q[1];
rz(-0.15730102) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0387715) q[0];
sx q[0];
rz(-0.73909315) q[0];
sx q[0];
rz(-3.0022013) q[0];
rz(-pi) q[1];
rz(0.83432014) q[2];
sx q[2];
rz(-1.7322707) q[2];
sx q[2];
rz(1.8404567) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3992608) q[1];
sx q[1];
rz(-2.4943588) q[1];
sx q[1];
rz(1.8982235) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62742426) q[3];
sx q[3];
rz(-1.809123) q[3];
sx q[3];
rz(2.0605007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9592231) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(-2.3106993) q[2];
rz(-0.043878555) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(-0.20251814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4956932) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(-0.010852531) q[0];
rz(0.27663484) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(0.29327926) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0115259) q[0];
sx q[0];
rz(-1.5767326) q[0];
sx q[0];
rz(2.5991711) q[0];
rz(-pi) q[1];
rz(1.3515527) q[2];
sx q[2];
rz(-0.65537894) q[2];
sx q[2];
rz(3.0041681) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5756702) q[1];
sx q[1];
rz(-2.2568963) q[1];
sx q[1];
rz(2.5887606) q[1];
rz(-pi) q[2];
rz(2.6777612) q[3];
sx q[3];
rz(-1.2004735) q[3];
sx q[3];
rz(-2.2390389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.098112) q[2];
sx q[2];
rz(-0.094823368) q[2];
sx q[2];
rz(-3.0528255) q[2];
rz(-0.25012112) q[3];
sx q[3];
rz(-1.1238031) q[3];
sx q[3];
rz(0.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1373238) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(-2.0776757) q[0];
rz(2.6990199) q[1];
sx q[1];
rz(-1.4241709) q[1];
sx q[1];
rz(1.3508505) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3063072) q[0];
sx q[0];
rz(-1.4915691) q[0];
sx q[0];
rz(-1.6552734) q[0];
rz(0.66019085) q[2];
sx q[2];
rz(-1.635951) q[2];
sx q[2];
rz(0.24812631) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6752154) q[1];
sx q[1];
rz(-1.5870759) q[1];
sx q[1];
rz(-2.0051763) q[1];
rz(-pi) q[2];
rz(-3.0382022) q[3];
sx q[3];
rz(-1.7915465) q[3];
sx q[3];
rz(3.0605928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.110048) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(2.6943977) q[2];
rz(-1.3859008) q[3];
sx q[3];
rz(-0.30062506) q[3];
sx q[3];
rz(-0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8266325) q[0];
sx q[0];
rz(-1.5216014) q[0];
sx q[0];
rz(0.78053027) q[0];
rz(-2.7138117) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(-0.16960493) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1113838) q[0];
sx q[0];
rz(-1.9281403) q[0];
sx q[0];
rz(1.0650728) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5461966) q[2];
sx q[2];
rz(-1.9925756) q[2];
sx q[2];
rz(-1.3299024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6360213) q[1];
sx q[1];
rz(-2.1850056) q[1];
sx q[1];
rz(-3.0843656) q[1];
rz(-1.9785089) q[3];
sx q[3];
rz(-2.793503) q[3];
sx q[3];
rz(-0.40504328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2910989) q[2];
sx q[2];
rz(-1.9591745) q[2];
sx q[2];
rz(-0.69236857) q[2];
rz(-2.678357) q[3];
sx q[3];
rz(-1.4866456) q[3];
sx q[3];
rz(2.0326116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8511843) q[0];
sx q[0];
rz(-1.52607) q[0];
sx q[0];
rz(-2.4551256) q[0];
rz(-2.6295173) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(1.2693263) q[2];
sx q[2];
rz(-1.8282679) q[2];
sx q[2];
rz(-1.4399583) q[2];
rz(2.0797707) q[3];
sx q[3];
rz(-2.90425) q[3];
sx q[3];
rz(-1.4958285) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];