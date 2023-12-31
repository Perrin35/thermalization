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
rz(0.16103345) q[0];
rz(1.610202) q[1];
sx q[1];
rz(-0.47709563) q[1];
sx q[1];
rz(0.49638003) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.046145) q[0];
sx q[0];
rz(-0.60337043) q[0];
sx q[0];
rz(-1.3674111) q[0];
rz(-2.8018087) q[2];
sx q[2];
rz(-1.7614363) q[2];
sx q[2];
rz(-1.9356188) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1094184) q[1];
sx q[1];
rz(-0.86583455) q[1];
sx q[1];
rz(1.6900307) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7339488) q[3];
sx q[3];
rz(-1.4789494) q[3];
sx q[3];
rz(-0.44029217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6266142) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(0.28764763) q[2];
rz(1.7488165) q[3];
sx q[3];
rz(-2.1578622) q[3];
sx q[3];
rz(-1.1028517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87067938) q[0];
sx q[0];
rz(-0.7551071) q[0];
sx q[0];
rz(-2.5640008) q[0];
rz(-1.6060991) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(1.5637406) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1801115) q[0];
sx q[0];
rz(-2.8948088) q[0];
sx q[0];
rz(-0.47829511) q[0];
rz(-pi) q[1];
rz(-0.30303843) q[2];
sx q[2];
rz(-1.8327456) q[2];
sx q[2];
rz(2.9532202) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8375081) q[1];
sx q[1];
rz(-1.6581074) q[1];
sx q[1];
rz(1.1156032) q[1];
rz(0.79264499) q[3];
sx q[3];
rz(-2.9229197) q[3];
sx q[3];
rz(1.6817026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0343798) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(2.0593026) q[2];
rz(1.1632495) q[3];
sx q[3];
rz(-1.2000368) q[3];
sx q[3];
rz(-0.64003402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0619693) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(2.9586155) q[0];
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
rz(-0.60488932) q[0];
sx q[0];
rz(-2.0920144) q[0];
sx q[0];
rz(-2.6813566) q[0];
rz(-pi) q[1];
x q[1];
rz(1.794533) q[2];
sx q[2];
rz(-1.1615331) q[2];
sx q[2];
rz(0.70659107) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.548) q[1];
sx q[1];
rz(-1.2263733) q[1];
sx q[1];
rz(-2.6542632) q[1];
rz(0.010300962) q[3];
sx q[3];
rz(-0.28763887) q[3];
sx q[3];
rz(1.079263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1220876) q[2];
sx q[2];
rz(-1.1818037) q[2];
sx q[2];
rz(2.2369177) q[2];
rz(2.4217862) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(-2.3864746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(1.6782137) q[0];
sx q[0];
rz(-0.97020522) q[0];
sx q[0];
rz(2.4257207) q[0];
rz(-0.32575682) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(2.148927) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7870165) q[0];
sx q[0];
rz(-1.0542608) q[0];
sx q[0];
rz(1.2660962) q[0];
x q[1];
rz(0.86187141) q[2];
sx q[2];
rz(-2.3700691) q[2];
sx q[2];
rz(2.9574403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.042295481) q[1];
sx q[1];
rz(-0.18810939) q[1];
sx q[1];
rz(0.52109615) q[1];
x q[2];
rz(-2.1487001) q[3];
sx q[3];
rz(-1.9979457) q[3];
sx q[3];
rz(-0.072673365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0585534) q[2];
sx q[2];
rz(-0.71648592) q[2];
sx q[2];
rz(-1.7129718) q[2];
rz(-0.84609091) q[3];
sx q[3];
rz(-1.5139791) q[3];
sx q[3];
rz(1.9184453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5457526) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(0.76675057) q[0];
rz(1.9013566) q[1];
sx q[1];
rz(-0.84754544) q[1];
sx q[1];
rz(-0.3516745) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7809217) q[0];
sx q[0];
rz(-2.43649) q[0];
sx q[0];
rz(2.6893925) q[0];
rz(-pi) q[1];
rz(2.211116) q[2];
sx q[2];
rz(-0.74379197) q[2];
sx q[2];
rz(3.004068) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.60492086) q[1];
sx q[1];
rz(-2.7465638) q[1];
sx q[1];
rz(0.037996304) q[1];
rz(-pi) q[2];
rz(1.0655754) q[3];
sx q[3];
rz(-1.8540283) q[3];
sx q[3];
rz(-3.1056044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9436283) q[2];
sx q[2];
rz(-1.6683656) q[2];
sx q[2];
rz(-2.6082805) q[2];
rz(2.7126281) q[3];
sx q[3];
rz(-0.52754378) q[3];
sx q[3];
rz(-2.0146446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-2.922373) q[1];
sx q[1];
rz(2.0419962) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0352286) q[0];
sx q[0];
rz(-1.117825) q[0];
sx q[0];
rz(-0.41082541) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9526688) q[2];
sx q[2];
rz(-1.0475698) q[2];
sx q[2];
rz(-0.97666937) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1814327) q[1];
sx q[1];
rz(-2.2391041) q[1];
sx q[1];
rz(-1.3431576) q[1];
rz(-1.908329) q[3];
sx q[3];
rz(-0.81157717) q[3];
sx q[3];
rz(2.8525762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.6301443) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(2.8586094) q[2];
rz(2.4354368) q[3];
sx q[3];
rz(-1.5420087) q[3];
sx q[3];
rz(-0.63505665) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0866078) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(-1.6116066) q[0];
rz(0.64487547) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(-2.9842916) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8511714) q[0];
sx q[0];
rz(-0.84050814) q[0];
sx q[0];
rz(1.4448326) q[0];
x q[1];
rz(-0.83432014) q[2];
sx q[2];
rz(-1.4093219) q[2];
sx q[2];
rz(1.8404567) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7423319) q[1];
sx q[1];
rz(-2.4943588) q[1];
sx q[1];
rz(-1.2433692) q[1];
x q[2];
rz(0.62742426) q[3];
sx q[3];
rz(-1.3324696) q[3];
sx q[3];
rz(2.0605007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1823696) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(0.83089337) q[2];
rz(-0.043878555) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(-0.20251814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6458994) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(0.010852531) q[0];
rz(2.8649578) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(-0.29327926) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5921708) q[0];
sx q[0];
rz(-0.54245078) q[0];
sx q[0];
rz(0.011499238) q[0];
x q[1];
rz(2.2145055) q[2];
sx q[2];
rz(-1.4378528) q[2];
sx q[2];
rz(1.5333652) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.20626942) q[1];
sx q[1];
rz(-2.2895797) q[1];
sx q[1];
rz(-2.1410336) q[1];
x q[2];
rz(1.9803489) q[3];
sx q[3];
rz(-1.1405986) q[3];
sx q[3];
rz(-0.48914117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.098112) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(3.0528255) q[2];
rz(-0.25012112) q[3];
sx q[3];
rz(-2.0177896) q[3];
sx q[3];
rz(-0.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0042689) q[0];
sx q[0];
rz(-0.22432888) q[0];
sx q[0];
rz(1.0639169) q[0];
rz(-2.6990199) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(1.3508505) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3063072) q[0];
sx q[0];
rz(-1.4915691) q[0];
sx q[0];
rz(-1.6552734) q[0];
x q[1];
rz(2.4814018) q[2];
sx q[2];
rz(-1.635951) q[2];
sx q[2];
rz(2.8934663) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6752154) q[1];
sx q[1];
rz(-1.5545168) q[1];
sx q[1];
rz(2.0051763) q[1];
rz(-pi) q[2];
rz(-2.001858) q[3];
sx q[3];
rz(-2.8981879) q[3];
sx q[3];
rz(-2.7800625) q[3];
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
rz(0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31496012) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(-2.3610624) q[0];
rz(-2.7138117) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(2.9719877) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7923332) q[0];
sx q[0];
rz(-1.0997286) q[0];
sx q[0];
rz(-0.40339289) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7197003) q[2];
sx q[2];
rz(-1.5932398) q[2];
sx q[2];
rz(2.9107712) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5055713) q[1];
sx q[1];
rz(-0.95658703) q[1];
sx q[1];
rz(0.057227055) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8923558) q[3];
sx q[3];
rz(-1.7064629) q[3];
sx q[3];
rz(-0.78007573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85049373) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(-0.69236857) q[2];
rz(2.678357) q[3];
sx q[3];
rz(-1.4866456) q[3];
sx q[3];
rz(1.108981) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(1.7789755) q[3];
sx q[3];
rz(-1.6856185) q[3];
sx q[3];
rz(0.57193397) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
