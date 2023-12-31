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
rz(2.0523235) q[0];
sx q[0];
rz(9.2637445) q[0];
rz(1.610202) q[1];
sx q[1];
rz(-0.47709563) q[1];
sx q[1];
rz(-2.6452126) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.046145) q[0];
sx q[0];
rz(-2.5382222) q[0];
sx q[0];
rz(-1.3674111) q[0];
rz(-pi) q[1];
rz(-0.52486323) q[2];
sx q[2];
rz(-0.3877936) q[2];
sx q[2];
rz(2.2847069) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1094184) q[1];
sx q[1];
rz(-0.86583455) q[1];
sx q[1];
rz(1.6900307) q[1];
rz(0.40764383) q[3];
sx q[3];
rz(-1.6626433) q[3];
sx q[3];
rz(-2.7013005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.51497841) q[2];
sx q[2];
rz(-0.86003059) q[2];
sx q[2];
rz(2.853945) q[2];
rz(1.3927762) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(2.0387409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2709133) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(-2.5640008) q[0];
rz(1.6060991) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(1.577852) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0664415) q[0];
sx q[0];
rz(-1.6834714) q[0];
sx q[0];
rz(0.22002797) q[0];
x q[1];
rz(-0.73194506) q[2];
sx q[2];
rz(-2.7436896) q[2];
sx q[2];
rz(2.0741472) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.3040846) q[1];
sx q[1];
rz(-1.4834852) q[1];
sx q[1];
rz(-1.1156032) q[1];
x q[2];
rz(-2.3489477) q[3];
sx q[3];
rz(-2.9229197) q[3];
sx q[3];
rz(1.6817026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0343798) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(1.0822901) q[2];
rz(1.1632495) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(0.64003402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0796233) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(2.9586155) q[0];
rz(-3.0984763) q[1];
sx q[1];
rz(-0.95298195) q[1];
sx q[1];
rz(-1.144369) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9632918) q[0];
sx q[0];
rz(-2.460647) q[0];
sx q[0];
rz(0.91239022) q[0];
x q[1];
rz(1.794533) q[2];
sx q[2];
rz(-1.9800595) q[2];
sx q[2];
rz(2.4350016) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.548) q[1];
sx q[1];
rz(-1.9152194) q[1];
sx q[1];
rz(0.48732948) q[1];
rz(-pi) q[2];
x q[2];
rz(0.010300962) q[3];
sx q[3];
rz(-2.8539538) q[3];
sx q[3];
rz(2.0623296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.019505067) q[2];
sx q[2];
rz(-1.1818037) q[2];
sx q[2];
rz(0.90467492) q[2];
rz(-2.4217862) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(-0.75511801) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6782137) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(-2.4257207) q[0];
rz(2.8158358) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(0.99266565) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216114) q[0];
sx q[0];
rz(-0.5926168) q[0];
sx q[0];
rz(0.48595925) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2797212) q[2];
sx q[2];
rz(-0.77152354) q[2];
sx q[2];
rz(2.9574403) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.48651931) q[1];
sx q[1];
rz(-1.4078948) q[1];
sx q[1];
rz(-1.4763114) q[1];
x q[2];
rz(-0.99289258) q[3];
sx q[3];
rz(-1.9979457) q[3];
sx q[3];
rz(-3.0689193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0830393) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(-1.7129718) q[2];
rz(-2.2955017) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(-1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59584004) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(-2.3748421) q[0];
rz(-1.2402361) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(-2.7899182) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2131166) q[0];
sx q[0];
rz(-0.94841829) q[0];
sx q[0];
rz(-1.2147796) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.211116) q[2];
sx q[2];
rz(-0.74379197) q[2];
sx q[2];
rz(-3.004068) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1406447) q[1];
sx q[1];
rz(-1.5854156) q[1];
sx q[1];
rz(-0.39477243) q[1];
x q[2];
rz(0.3211) q[3];
sx q[3];
rz(-2.0541111) q[3];
sx q[3];
rz(1.4534284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.19796431) q[2];
sx q[2];
rz(-1.4732271) q[2];
sx q[2];
rz(-0.53331214) q[2];
rz(-0.42896459) q[3];
sx q[3];
rz(-2.6140489) q[3];
sx q[3];
rz(-1.126948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11739843) q[0];
sx q[0];
rz(-2.0485931) q[0];
sx q[0];
rz(2.5999516) q[0];
rz(-0.60846865) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(2.0419962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8894419) q[0];
sx q[0];
rz(-0.60177207) q[0];
sx q[0];
rz(0.88366951) q[0];
rz(-pi) q[1];
rz(-1.2560558) q[2];
sx q[2];
rz(-0.55329269) q[2];
sx q[2];
rz(1.3421343) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.82367831) q[1];
sx q[1];
rz(-0.70033973) q[1];
sx q[1];
rz(-2.8631696) q[1];
rz(1.908329) q[3];
sx q[3];
rz(-0.81157717) q[3];
sx q[3];
rz(-2.8525762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5114484) q[2];
sx q[2];
rz(-1.6161329) q[2];
sx q[2];
rz(0.28298322) q[2];
rz(-2.4354368) q[3];
sx q[3];
rz(-1.5420087) q[3];
sx q[3];
rz(0.63505665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-0.054984897) q[0];
sx q[0];
rz(-0.98751706) q[0];
sx q[0];
rz(1.5299861) q[0];
rz(2.4967172) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(2.9842916) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29042127) q[0];
sx q[0];
rz(-0.84050814) q[0];
sx q[0];
rz(-1.4448326) q[0];
x q[1];
rz(-0.83432014) q[2];
sx q[2];
rz(-1.4093219) q[2];
sx q[2];
rz(-1.3011359) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7423319) q[1];
sx q[1];
rz(-2.4943588) q[1];
sx q[1];
rz(1.2433692) q[1];
x q[2];
rz(-1.2792475) q[3];
sx q[3];
rz(-2.1778717) q[3];
sx q[3];
rz(0.32015043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9592231) q[2];
sx q[2];
rz(-0.59252858) q[2];
sx q[2];
rz(2.3106993) q[2];
rz(3.0977141) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(2.9390745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6458994) q[0];
sx q[0];
rz(-0.93911397) q[0];
sx q[0];
rz(3.1307401) q[0];
rz(-0.27663484) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(2.8483134) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56284833) q[0];
sx q[0];
rz(-1.0283854) q[0];
sx q[0];
rz(1.5777274) q[0];
x q[1];
rz(0.92708712) q[2];
sx q[2];
rz(-1.4378528) q[2];
sx q[2];
rz(1.6082275) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3775023) q[1];
sx q[1];
rz(-1.9891771) q[1];
sx q[1];
rz(2.3368895) q[1];
rz(-0.71473748) q[3];
sx q[3];
rz(-0.5849896) q[3];
sx q[3];
rz(-1.294567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.098112) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(3.0528255) q[2];
rz(2.8914715) q[3];
sx q[3];
rz(-2.0177896) q[3];
sx q[3];
rz(2.5626101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0042689) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(1.0639169) q[0];
rz(-0.44257277) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(-1.3508505) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0161184) q[0];
sx q[0];
rz(-0.11575143) q[0];
sx q[0];
rz(2.3257757) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4883792) q[2];
sx q[2];
rz(-2.2293408) q[2];
sx q[2];
rz(1.272162) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0722479) q[1];
sx q[1];
rz(-2.7069271) q[1];
sx q[1];
rz(1.6094631) q[1];
rz(-1.3488995) q[3];
sx q[3];
rz(-1.4699234) q[3];
sx q[3];
rz(-1.6290806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0315447) q[2];
sx q[2];
rz(-1.2908547) q[2];
sx q[2];
rz(2.6943977) q[2];
rz(1.7556919) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(-2.6269004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31496012) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(2.3610624) q[0];
rz(-2.7138117) q[1];
sx q[1];
rz(-2.784274) q[1];
sx q[1];
rz(-2.9719877) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34925941) q[0];
sx q[0];
rz(-1.0997286) q[0];
sx q[0];
rz(-0.40339289) q[0];
rz(-pi) q[1];
rz(-1.5461966) q[2];
sx q[2];
rz(-1.149017) q[2];
sx q[2];
rz(1.8116902) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6360213) q[1];
sx q[1];
rz(-2.1850056) q[1];
sx q[1];
rz(-3.0843656) q[1];
rz(-pi) q[2];
rz(-1.9785089) q[3];
sx q[3];
rz(-2.793503) q[3];
sx q[3];
rz(2.7365494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2910989) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(0.69236857) q[2];
rz(0.46323562) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(-2.0326116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2904084) q[0];
sx q[0];
rz(-1.6155227) q[0];
sx q[0];
rz(0.68646705) q[0];
rz(-0.51207536) q[1];
sx q[1];
rz(-2.8072186) q[1];
sx q[1];
rz(1.0288815) q[1];
rz(-2.2962773) q[2];
sx q[2];
rz(-2.7477063) q[2];
sx q[2];
rz(0.817) q[2];
rz(3.0242596) q[3];
sx q[3];
rz(-1.7775848) q[3];
sx q[3];
rz(2.1669273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
