OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.53192294) q[0];
sx q[0];
rz(1.6165531) q[0];
sx q[0];
rz(11.308148) q[0];
rz(-1.582107) q[1];
sx q[1];
rz(-2.6546302) q[1];
sx q[1];
rz(2.7064986) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5814911) q[0];
sx q[0];
rz(-1.035443) q[0];
sx q[0];
rz(2.125432) q[0];
x q[1];
rz(1.9149532) q[2];
sx q[2];
rz(-1.0593057) q[2];
sx q[2];
rz(2.34735) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62578979) q[1];
sx q[1];
rz(-0.28072673) q[1];
sx q[1];
rz(1.8516385) q[1];
rz(-pi) q[2];
rz(-2.4354127) q[3];
sx q[3];
rz(-0.67863388) q[3];
sx q[3];
rz(0.61509252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2121928) q[2];
sx q[2];
rz(-1.2779028) q[2];
sx q[2];
rz(-1.0431935) q[2];
rz(0.4010525) q[3];
sx q[3];
rz(-1.6845104) q[3];
sx q[3];
rz(-0.57798398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9143518) q[0];
sx q[0];
rz(-2.9965017) q[0];
sx q[0];
rz(-2.5703854) q[0];
rz(-0.97194833) q[1];
sx q[1];
rz(-2.4010039) q[1];
sx q[1];
rz(-2.0170225) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49405022) q[0];
sx q[0];
rz(-2.4801773) q[0];
sx q[0];
rz(-1.1986092) q[0];
rz(-2.6827181) q[2];
sx q[2];
rz(-2.4374408) q[2];
sx q[2];
rz(1.6353232) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.827885) q[1];
sx q[1];
rz(-1.4769625) q[1];
sx q[1];
rz(-0.83637823) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5519106) q[3];
sx q[3];
rz(-1.369993) q[3];
sx q[3];
rz(2.7310179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9198415) q[2];
sx q[2];
rz(-1.5849179) q[2];
sx q[2];
rz(2.705503) q[2];
rz(-0.74887216) q[3];
sx q[3];
rz(-0.80513969) q[3];
sx q[3];
rz(2.3124636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1083168) q[0];
sx q[0];
rz(-1.746614) q[0];
sx q[0];
rz(-3.0535611) q[0];
rz(2.2875359) q[1];
sx q[1];
rz(-0.66103926) q[1];
sx q[1];
rz(1.0850272) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4861761) q[0];
sx q[0];
rz(-1.2987505) q[0];
sx q[0];
rz(-0.06975091) q[0];
rz(2.9879007) q[2];
sx q[2];
rz(-1.7847381) q[2];
sx q[2];
rz(-0.74665507) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.67133707) q[1];
sx q[1];
rz(-0.47154676) q[1];
sx q[1];
rz(1.9370473) q[1];
rz(-pi) q[2];
rz(0.17739968) q[3];
sx q[3];
rz(-2.7466603) q[3];
sx q[3];
rz(1.653914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.4714841) q[2];
sx q[2];
rz(-1.4205616) q[2];
sx q[2];
rz(-2.3954083) q[2];
rz(-0.21607312) q[3];
sx q[3];
rz(-1.2553299) q[3];
sx q[3];
rz(0.20572534) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7699319) q[0];
sx q[0];
rz(-0.07337229) q[0];
sx q[0];
rz(1.7727456) q[0];
rz(-0.61028496) q[1];
sx q[1];
rz(-1.7758324) q[1];
sx q[1];
rz(2.761421) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1606522) q[0];
sx q[0];
rz(-1.9987172) q[0];
sx q[0];
rz(-1.9141657) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39003464) q[2];
sx q[2];
rz(-0.48214285) q[2];
sx q[2];
rz(0.48393656) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4277162) q[1];
sx q[1];
rz(-2.6170934) q[1];
sx q[1];
rz(-1.8668411) q[1];
rz(-pi) q[2];
rz(0.065297619) q[3];
sx q[3];
rz(-1.9717798) q[3];
sx q[3];
rz(0.15006615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.88468203) q[2];
sx q[2];
rz(-0.94748777) q[2];
sx q[2];
rz(-1.0379637) q[2];
rz(2.7773618) q[3];
sx q[3];
rz(-2.3528152) q[3];
sx q[3];
rz(-2.4558333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6895741) q[0];
sx q[0];
rz(-1.4301393) q[0];
sx q[0];
rz(0.83056393) q[0];
rz(0.05016249) q[1];
sx q[1];
rz(-3.0144033) q[1];
sx q[1];
rz(-0.62279472) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2723214) q[0];
sx q[0];
rz(-1.4928482) q[0];
sx q[0];
rz(2.9236631) q[0];
x q[1];
rz(-0.41983126) q[2];
sx q[2];
rz(-11/(2*pi)) q[2];
sx q[2];
rz(-0.77984389) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2166136) q[1];
sx q[1];
rz(-2.2053524) q[1];
sx q[1];
rz(-2.6960009) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5708617) q[3];
sx q[3];
rz(-1.5932461) q[3];
sx q[3];
rz(2.8714436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.043776) q[2];
sx q[2];
rz(-2.706683) q[2];
sx q[2];
rz(-0.80769509) q[2];
rz(1.5378599) q[3];
sx q[3];
rz(-1.6350428) q[3];
sx q[3];
rz(2.9174771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1767126) q[0];
sx q[0];
rz(-1.3310615) q[0];
sx q[0];
rz(1.4759195) q[0];
rz(1.2337947) q[1];
sx q[1];
rz(-1.7314311) q[1];
sx q[1];
rz(2.9880611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3672402) q[0];
sx q[0];
rz(-1.4611547) q[0];
sx q[0];
rz(1.5745274) q[0];
rz(2.0664882) q[2];
sx q[2];
rz(-1.1678936) q[2];
sx q[2];
rz(-2.7893381) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1099838) q[1];
sx q[1];
rz(-0.69797198) q[1];
sx q[1];
rz(-2.4137761) q[1];
rz(2.5028822) q[3];
sx q[3];
rz(-1.6473466) q[3];
sx q[3];
rz(-2.9552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.70225707) q[2];
sx q[2];
rz(-2.3372237) q[2];
sx q[2];
rz(-0.55629936) q[2];
rz(3.124681) q[3];
sx q[3];
rz(-0.97512475) q[3];
sx q[3];
rz(-2.4771966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9282114) q[0];
sx q[0];
rz(-0.072755486) q[0];
sx q[0];
rz(0.67759204) q[0];
rz(-1.821359) q[1];
sx q[1];
rz(-1.0878539) q[1];
sx q[1];
rz(-2.5755612) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6460577) q[0];
sx q[0];
rz(-1.2355775) q[0];
sx q[0];
rz(2.6779102) q[0];
x q[1];
rz(1.2957057) q[2];
sx q[2];
rz(-0.054271532) q[2];
sx q[2];
rz(-1.4473297) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.89443242) q[1];
sx q[1];
rz(-1.4196779) q[1];
sx q[1];
rz(-2.2770348) q[1];
rz(2.651068) q[3];
sx q[3];
rz(-1.5340866) q[3];
sx q[3];
rz(1.6877268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2322106) q[2];
sx q[2];
rz(-2.7558694) q[2];
sx q[2];
rz(-2.5717226) q[2];
rz(-1.0424403) q[3];
sx q[3];
rz(-1.9614377) q[3];
sx q[3];
rz(2.1377783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7637699) q[0];
sx q[0];
rz(-0.42709392) q[0];
sx q[0];
rz(-1.0231934) q[0];
rz(-2.2238253) q[1];
sx q[1];
rz(-2.387391) q[1];
sx q[1];
rz(0.86407152) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9167332) q[0];
sx q[0];
rz(-1.8141439) q[0];
sx q[0];
rz(-0.12638522) q[0];
rz(-pi) q[1];
rz(2.7644453) q[2];
sx q[2];
rz(-2.0066064) q[2];
sx q[2];
rz(-0.47919861) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.041124972) q[1];
sx q[1];
rz(-2.0534614) q[1];
sx q[1];
rz(-2.3339999) q[1];
x q[2];
rz(-0.72743639) q[3];
sx q[3];
rz(-1.2890649) q[3];
sx q[3];
rz(1.7121332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0495905) q[2];
sx q[2];
rz(-1.8747753) q[2];
sx q[2];
rz(-0.76466307) q[2];
rz(-2.2406254) q[3];
sx q[3];
rz(-0.67470208) q[3];
sx q[3];
rz(-1.0990134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6242999) q[0];
sx q[0];
rz(-2.7462672) q[0];
sx q[0];
rz(0.98315352) q[0];
rz(-0.82955018) q[1];
sx q[1];
rz(-1.364565) q[1];
sx q[1];
rz(0.57873908) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0522636) q[0];
sx q[0];
rz(-2.7476285) q[0];
sx q[0];
rz(1.5376505) q[0];
rz(-pi) q[1];
rz(1.0320682) q[2];
sx q[2];
rz(-1.0916289) q[2];
sx q[2];
rz(0.62512809) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1956801) q[1];
sx q[1];
rz(-2.0766313) q[1];
sx q[1];
rz(-2.6786184) q[1];
rz(0.10161256) q[3];
sx q[3];
rz(-0.34694537) q[3];
sx q[3];
rz(-0.24747758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4911554) q[2];
sx q[2];
rz(-1.8530242) q[2];
sx q[2];
rz(3.0299661) q[2];
rz(-0.95587436) q[3];
sx q[3];
rz(-0.99146944) q[3];
sx q[3];
rz(-1.2607695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2535506) q[0];
sx q[0];
rz(-1.0567867) q[0];
sx q[0];
rz(0.014130935) q[0];
rz(-1.1125394) q[1];
sx q[1];
rz(-2.2281149) q[1];
sx q[1];
rz(1.5412451) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80373549) q[0];
sx q[0];
rz(-1.4232262) q[0];
sx q[0];
rz(-1.4413367) q[0];
rz(-3.0156323) q[2];
sx q[2];
rz(-2.2114814) q[2];
sx q[2];
rz(1.0116742) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4804876) q[1];
sx q[1];
rz(-0.36927642) q[1];
sx q[1];
rz(-2.9758456) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16997204) q[3];
sx q[3];
rz(-1.7331496) q[3];
sx q[3];
rz(1.270783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17497002) q[2];
sx q[2];
rz(-1.2302159) q[2];
sx q[2];
rz(2.5959385) q[2];
rz(-3.0555861) q[3];
sx q[3];
rz(-1.6274446) q[3];
sx q[3];
rz(-2.8470794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1957112) q[0];
sx q[0];
rz(-2.0989037) q[0];
sx q[0];
rz(1.209191) q[0];
rz(-0.60278268) q[1];
sx q[1];
rz(-0.61059112) q[1];
sx q[1];
rz(0.75377725) q[1];
rz(-2.9309737) q[2];
sx q[2];
rz(-0.89952614) q[2];
sx q[2];
rz(2.6256403) q[2];
rz(0.62282563) q[3];
sx q[3];
rz(-1.4824184) q[3];
sx q[3];
rz(-1.3893954) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
