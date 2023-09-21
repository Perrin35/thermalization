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
rz(2.664497) q[1];
sx q[1];
rz(8.9283979) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8007322) q[0];
sx q[0];
rz(-2.1600318) q[0];
sx q[0];
rz(-0.13829921) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33978396) q[2];
sx q[2];
rz(-1.7614363) q[2];
sx q[2];
rz(-1.9356188) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5254933) q[1];
sx q[1];
rz(-1.4800737) q[1];
sx q[1];
rz(-0.70848042) q[1];
rz(-pi) q[2];
rz(-0.22827893) q[3];
sx q[3];
rz(-2.7242959) q[3];
sx q[3];
rz(-1.8018064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51497841) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(-0.28764763) q[2];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87067938) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(-2.5640008) q[0];
rz(-1.6060991) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(1.5637406) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47050414) q[0];
sx q[0];
rz(-1.3521863) q[0];
sx q[0];
rz(1.455362) q[0];
x q[1];
rz(0.30303843) q[2];
sx q[2];
rz(-1.3088471) q[2];
sx q[2];
rz(-0.18837243) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2240552) q[1];
sx q[1];
rz(-1.1174669) q[1];
sx q[1];
rz(0.097150306) q[1];
rz(0.15474774) q[3];
sx q[3];
rz(-1.7259211) q[3];
sx q[3];
rz(-2.250092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0343798) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(-1.0822901) q[2];
rz(-1.9783431) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(0.64003402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0796233) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(2.9586155) q[0];
rz(3.0984763) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(1.9972237) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60488932) q[0];
sx q[0];
rz(-2.0920144) q[0];
sx q[0];
rz(0.46023603) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3470596) q[2];
sx q[2];
rz(-1.9800595) q[2];
sx q[2];
rz(-2.4350016) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3414587) q[1];
sx q[1];
rz(-1.1143436) q[1];
sx q[1];
rz(1.1851428) q[1];
rz(0.010300962) q[3];
sx q[3];
rz(-0.28763887) q[3];
sx q[3];
rz(-2.0623296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1220876) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(-2.2369177) q[2];
rz(-0.71980643) q[3];
sx q[3];
rz(-2.7147839) q[3];
sx q[3];
rz(-0.75511801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4633789) q[0];
sx q[0];
rz(-0.97020522) q[0];
sx q[0];
rz(2.4257207) q[0];
rz(-0.32575682) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(2.148927) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2199812) q[0];
sx q[0];
rz(-2.5489759) q[0];
sx q[0];
rz(-0.48595925) q[0];
x q[1];
rz(2.2797212) q[2];
sx q[2];
rz(-0.77152354) q[2];
sx q[2];
rz(-0.18415235) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0996453) q[1];
sx q[1];
rz(-1.6640267) q[1];
sx q[1];
rz(-0.16361841) q[1];
rz(-pi) q[2];
rz(2.2654815) q[3];
sx q[3];
rz(-2.4377341) q[3];
sx q[3];
rz(1.0775523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0585534) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(-1.4286208) q[2];
rz(-2.2955017) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(-1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59584004) q[0];
sx q[0];
rz(-1.4056453) q[0];
sx q[0];
rz(0.76675057) q[0];
rz(1.9013566) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(-2.7899182) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2857916) q[0];
sx q[0];
rz(-1.283678) q[0];
sx q[0];
rz(2.48824) q[0];
rz(2.211116) q[2];
sx q[2];
rz(-2.3978007) q[2];
sx q[2];
rz(-3.004068) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5778351) q[1];
sx q[1];
rz(-1.1760684) q[1];
sx q[1];
rz(1.5866336) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1122123) q[3];
sx q[3];
rz(-2.5684528) q[3];
sx q[3];
rz(2.074632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9436283) q[2];
sx q[2];
rz(-1.4732271) q[2];
sx q[2];
rz(2.6082805) q[2];
rz(2.7126281) q[3];
sx q[3];
rz(-0.52754378) q[3];
sx q[3];
rz(-2.0146446) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0241942) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(-2.5999516) q[0];
rz(-0.60846865) q[1];
sx q[1];
rz(-2.922373) q[1];
sx q[1];
rz(1.0995964) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.106364) q[0];
sx q[0];
rz(-2.0237676) q[0];
sx q[0];
rz(2.7307672) q[0];
x q[1];
rz(2.1018283) q[2];
sx q[2];
rz(-1.7341988) q[2];
sx q[2];
rz(0.49887564) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9601599) q[1];
sx q[1];
rz(-2.2391041) q[1];
sx q[1];
rz(-1.798435) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33575616) q[3];
sx q[3];
rz(-2.3246375) q[3];
sx q[3];
rz(-0.18248617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.6301443) q[2];
sx q[2];
rz(-1.6161329) q[2];
sx q[2];
rz(-2.8586094) q[2];
rz(0.7061559) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(-0.63505665) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054984897) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(1.6116066) q[0];
rz(-0.64487547) q[1];
sx q[1];
rz(-2.6480643) q[1];
sx q[1];
rz(-2.9842916) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.776942) q[0];
sx q[0];
rz(-1.4770664) q[0];
sx q[0];
rz(-0.73424299) q[0];
rz(-pi) q[1];
rz(-2.925161) q[2];
sx q[2];
rz(-2.2955403) q[2];
sx q[2];
rz(-3.016678) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0485222) q[1];
sx q[1];
rz(-1.375636) q[1];
sx q[1];
rz(-0.94961571) q[1];
x q[2];
rz(0.39237202) q[3];
sx q[3];
rz(-0.6654159) q[3];
sx q[3];
rz(0.80442807) q[3];
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
rz(2.9390745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.4956932) q[0];
sx q[0];
rz(-0.93911397) q[0];
sx q[0];
rz(-0.010852531) q[0];
rz(2.8649578) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(2.8483134) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56284833) q[0];
sx q[0];
rz(-2.1132073) q[0];
sx q[0];
rz(-1.5777274) q[0];
x q[1];
rz(2.2145055) q[2];
sx q[2];
rz(-1.4378528) q[2];
sx q[2];
rz(1.5333652) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.5659225) q[1];
sx q[1];
rz(-2.2568963) q[1];
sx q[1];
rz(-2.5887606) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46383143) q[3];
sx q[3];
rz(-1.9411191) q[3];
sx q[3];
rz(2.2390389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.04348065) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(0.088767178) q[2];
rz(-2.8914715) q[3];
sx q[3];
rz(-2.0177896) q[3];
sx q[3];
rz(-2.5626101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1373238) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(-2.0776757) q[0];
rz(-2.6990199) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(-1.7907422) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1254743) q[0];
sx q[0];
rz(-3.0258412) q[0];
sx q[0];
rz(2.3257757) q[0];
rz(-pi) q[1];
rz(-0.10599372) q[2];
sx q[2];
rz(-0.66291891) q[2];
sx q[2];
rz(1.7352599) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0693448) q[1];
sx q[1];
rz(-0.43466553) q[1];
sx q[1];
rz(-1.5321295) q[1];
rz(1.3488995) q[3];
sx q[3];
rz(-1.4699234) q[3];
sx q[3];
rz(1.6290806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0315447) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(-2.6943977) q[2];
rz(-1.7556919) q[3];
sx q[3];
rz(-0.30062506) q[3];
sx q[3];
rz(-2.6269004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31496012) q[0];
sx q[0];
rz(-1.5216014) q[0];
sx q[0];
rz(0.78053027) q[0];
rz(-0.42778095) q[1];
sx q[1];
rz(-2.784274) q[1];
sx q[1];
rz(-0.16960493) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34925941) q[0];
sx q[0];
rz(-2.0418641) q[0];
sx q[0];
rz(0.40339289) q[0];
x q[1];
rz(-0.42189235) q[2];
sx q[2];
rz(-1.5483529) q[2];
sx q[2];
rz(-2.9107712) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6360213) q[1];
sx q[1];
rz(-2.1850056) q[1];
sx q[1];
rz(0.057227055) q[1];
x q[2];
rz(1.8923558) q[3];
sx q[3];
rz(-1.4351298) q[3];
sx q[3];
rz(-2.3615169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.85049373) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(-0.69236857) q[2];
rz(0.46323562) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(-2.0326116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8511843) q[0];
sx q[0];
rz(-1.6155227) q[0];
sx q[0];
rz(0.68646705) q[0];
rz(-0.51207536) q[1];
sx q[1];
rz(-2.8072186) q[1];
sx q[1];
rz(1.0288815) q[1];
rz(-0.84531534) q[2];
sx q[2];
rz(-0.39388638) q[2];
sx q[2];
rz(-2.3245927) q[2];
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