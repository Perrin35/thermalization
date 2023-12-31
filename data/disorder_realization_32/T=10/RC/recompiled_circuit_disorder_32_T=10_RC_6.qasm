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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8344603) q[0];
sx q[0];
rz(-1.6856598) q[0];
sx q[0];
rz(-0.97712028) q[0];
rz(-1.7726937) q[2];
sx q[2];
rz(-1.9041833) q[2];
sx q[2];
rz(-0.29793973) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1094184) q[1];
sx q[1];
rz(-2.2757581) q[1];
sx q[1];
rz(-1.6900307) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4708038) q[3];
sx q[3];
rz(-1.9766207) q[3];
sx q[3];
rz(1.0909181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51497841) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(0.28764763) q[2];
rz(1.7488165) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(1.1028517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87067938) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(2.5640008) q[0];
rz(1.6060991) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(-1.5637406) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6710885) q[0];
sx q[0];
rz(-1.3521863) q[0];
sx q[0];
rz(1.455362) q[0];
x q[1];
rz(-0.30303843) q[2];
sx q[2];
rz(-1.3088471) q[2];
sx q[2];
rz(-2.9532202) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6985804) q[1];
sx q[1];
rz(-0.46291446) q[1];
sx q[1];
rz(-1.7673311) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9868449) q[3];
sx q[3];
rz(-1.7259211) q[3];
sx q[3];
rz(-2.250092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.10721283) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(2.0593026) q[2];
rz(-1.9783431) q[3];
sx q[3];
rz(-1.2000368) q[3];
sx q[3];
rz(2.5015586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0619693) q[0];
sx q[0];
rz(-0.29456961) q[0];
sx q[0];
rz(-0.18297718) q[0];
rz(0.043116365) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(1.144369) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5367033) q[0];
sx q[0];
rz(-2.0920144) q[0];
sx q[0];
rz(-0.46023603) q[0];
rz(-pi) q[1];
rz(-2.7230354) q[2];
sx q[2];
rz(-1.7757799) q[2];
sx q[2];
rz(2.1870854) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80013393) q[1];
sx q[1];
rz(-2.027249) q[1];
sx q[1];
rz(1.9564499) q[1];
rz(1.5677489) q[3];
sx q[3];
rz(-1.2831732) q[3];
sx q[3];
rz(-2.0730719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1220876) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(0.90467492) q[2];
rz(2.4217862) q[3];
sx q[3];
rz(-2.7147839) q[3];
sx q[3];
rz(2.3864746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4633789) q[0];
sx q[0];
rz(-0.97020522) q[0];
sx q[0];
rz(2.4257207) q[0];
rz(2.8158358) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(2.148927) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3545761) q[0];
sx q[0];
rz(-2.0873318) q[0];
sx q[0];
rz(-1.2660962) q[0];
rz(2.5771192) q[2];
sx q[2];
rz(-2.1285004) q[2];
sx q[2];
rz(-2.0828473) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.48651931) q[1];
sx q[1];
rz(-1.7336978) q[1];
sx q[1];
rz(-1.4763114) q[1];
rz(-pi) q[2];
rz(0.99289258) q[3];
sx q[3];
rz(-1.9979457) q[3];
sx q[3];
rz(3.0689193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0585534) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(1.7129718) q[2];
rz(-0.84609091) q[3];
sx q[3];
rz(-1.5139791) q[3];
sx q[3];
rz(1.9184453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59584004) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(-2.3748421) q[0];
rz(1.9013566) q[1];
sx q[1];
rz(-0.84754544) q[1];
sx q[1];
rz(2.7899182) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2857916) q[0];
sx q[0];
rz(-1.8579146) q[0];
sx q[0];
rz(-0.65335269) q[0];
rz(2.211116) q[2];
sx q[2];
rz(-2.3978007) q[2];
sx q[2];
rz(0.1375246) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5778351) q[1];
sx q[1];
rz(-1.9655242) q[1];
sx q[1];
rz(-1.5549591) q[1];
rz(-pi) q[2];
rz(-2.8204927) q[3];
sx q[3];
rz(-2.0541111) q[3];
sx q[3];
rz(1.4534284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.19796431) q[2];
sx q[2];
rz(-1.6683656) q[2];
sx q[2];
rz(2.6082805) q[2];
rz(-0.42896459) q[3];
sx q[3];
rz(-2.6140489) q[3];
sx q[3];
rz(2.0146446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11739843) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(2.5999516) q[0];
rz(-2.533124) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(1.0995964) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0352286) q[0];
sx q[0];
rz(-1.117825) q[0];
sx q[0];
rz(2.7307672) q[0];
x q[1];
rz(-1.0397644) q[2];
sx q[2];
rz(-1.7341988) q[2];
sx q[2];
rz(-2.642717) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.82367831) q[1];
sx q[1];
rz(-2.4412529) q[1];
sx q[1];
rz(-2.8631696) q[1];
rz(-1.2332637) q[3];
sx q[3];
rz(-2.3300155) q[3];
sx q[3];
rz(2.8525762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.6301443) q[2];
sx q[2];
rz(-1.6161329) q[2];
sx q[2];
rz(-2.8586094) q[2];
rz(2.4354368) q[3];
sx q[3];
rz(-1.5420087) q[3];
sx q[3];
rz(2.506536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0866078) q[0];
sx q[0];
rz(-0.98751706) q[0];
sx q[0];
rz(-1.6116066) q[0];
rz(2.4967172) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(2.9842916) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.776942) q[0];
sx q[0];
rz(-1.4770664) q[0];
sx q[0];
rz(-0.73424299) q[0];
x q[1];
rz(2.3072725) q[2];
sx q[2];
rz(-1.7322707) q[2];
sx q[2];
rz(1.3011359) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.09307043) q[1];
sx q[1];
rz(-1.7659566) q[1];
sx q[1];
rz(-0.94961571) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5141684) q[3];
sx q[3];
rz(-1.3324696) q[3];
sx q[3];
rz(1.081092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9592231) q[2];
sx q[2];
rz(-0.59252858) q[2];
sx q[2];
rz(-0.83089337) q[2];
rz(-3.0977141) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(-2.9390745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4956932) q[0];
sx q[0];
rz(-0.93911397) q[0];
sx q[0];
rz(0.010852531) q[0];
rz(-2.8649578) q[1];
sx q[1];
rz(-0.77181548) q[1];
sx q[1];
rz(-0.29327926) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1300668) q[0];
sx q[0];
rz(-1.5648601) q[0];
sx q[0];
rz(0.54242155) q[0];
rz(-pi) q[1];
rz(-0.16565928) q[2];
sx q[2];
rz(-0.93369166) q[2];
sx q[2];
rz(-0.13656244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.5659225) q[1];
sx q[1];
rz(-0.88469632) q[1];
sx q[1];
rz(-2.5887606) q[1];
x q[2];
rz(-2.4268552) q[3];
sx q[3];
rz(-0.5849896) q[3];
sx q[3];
rz(1.294567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.098112) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(-3.0528255) q[2];
rz(0.25012112) q[3];
sx q[3];
rz(-2.0177896) q[3];
sx q[3];
rz(0.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1373238) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(2.0776757) q[0];
rz(0.44257277) q[1];
sx q[1];
rz(-1.4241709) q[1];
sx q[1];
rz(-1.3508505) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1254743) q[0];
sx q[0];
rz(-0.11575143) q[0];
sx q[0];
rz(-2.3257757) q[0];
rz(-pi) q[1];
rz(2.4814018) q[2];
sx q[2];
rz(-1.5056416) q[2];
sx q[2];
rz(-2.8934663) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0296214) q[1];
sx q[1];
rz(-2.0051149) q[1];
sx q[1];
rz(-0.017945826) q[1];
x q[2];
rz(-0.10339046) q[3];
sx q[3];
rz(-1.7915465) q[3];
sx q[3];
rz(-3.0605928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.110048) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(2.6943977) q[2];
rz(1.3859008) q[3];
sx q[3];
rz(-0.30062506) q[3];
sx q[3];
rz(0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8266325) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(2.3610624) q[0];
rz(2.7138117) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(0.16960493) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0302089) q[0];
sx q[0];
rz(-1.9281403) q[0];
sx q[0];
rz(1.0650728) q[0];
rz(1.5461966) q[2];
sx q[2];
rz(-1.149017) q[2];
sx q[2];
rz(1.3299024) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.032220275) q[1];
sx q[1];
rz(-1.6175555) q[1];
sx q[1];
rz(-0.95581518) q[1];
rz(-1.2492368) q[3];
sx q[3];
rz(-1.7064629) q[3];
sx q[3];
rz(-0.78007573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.85049373) q[2];
sx q[2];
rz(-1.9591745) q[2];
sx q[2];
rz(2.4492241) q[2];
rz(-0.46323562) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(-1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2904084) q[0];
sx q[0];
rz(-1.52607) q[0];
sx q[0];
rz(-2.4551256) q[0];
rz(-2.6295173) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(-1.8722664) q[2];
sx q[2];
rz(-1.8282679) q[2];
sx q[2];
rz(-1.4399583) q[2];
rz(-3.0242596) q[3];
sx q[3];
rz(-1.3640079) q[3];
sx q[3];
rz(-0.97466536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
