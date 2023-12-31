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
rz(pi/2) q[2];
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
rz(2.8018087) q[2];
sx q[2];
rz(-1.7614363) q[2];
sx q[2];
rz(1.9356188) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1094184) q[1];
sx q[1];
rz(-2.2757581) q[1];
sx q[1];
rz(-1.4515619) q[1];
rz(0.22827893) q[3];
sx q[3];
rz(-0.41729673) q[3];
sx q[3];
rz(1.3397863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6266142) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(0.28764763) q[2];
rz(-1.3927762) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(1.1028517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2709133) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(0.57759181) q[0];
rz(1.6060991) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(1.577852) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0664415) q[0];
sx q[0];
rz(-1.4581212) q[0];
sx q[0];
rz(-0.22002797) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8385542) q[2];
sx q[2];
rz(-1.3088471) q[2];
sx q[2];
rz(-0.18837243) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2240552) q[1];
sx q[1];
rz(-2.0241258) q[1];
sx q[1];
rz(-3.0444423) q[1];
rz(-pi) q[2];
rz(0.79264499) q[3];
sx q[3];
rz(-0.21867293) q[3];
sx q[3];
rz(-1.6817026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.10721283) q[2];
sx q[2];
rz(-2.1652174) q[2];
sx q[2];
rz(2.0593026) q[2];
rz(-1.1632495) q[3];
sx q[3];
rz(-1.2000368) q[3];
sx q[3];
rz(0.64003402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0619693) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(-0.18297718) q[0];
rz(3.0984763) q[1];
sx q[1];
rz(-0.95298195) q[1];
sx q[1];
rz(1.144369) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9632918) q[0];
sx q[0];
rz(-2.460647) q[0];
sx q[0];
rz(2.2292024) q[0];
rz(-pi) q[1];
x q[1];
rz(1.794533) q[2];
sx q[2];
rz(-1.1615331) q[2];
sx q[2];
rz(0.70659107) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5444762) q[1];
sx q[1];
rz(-2.5529478) q[1];
sx q[1];
rz(-2.4878923) q[1];
rz(0.28762443) q[3];
sx q[3];
rz(-1.5678741) q[3];
sx q[3];
rz(-0.50141108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.019505067) q[2];
sx q[2];
rz(-1.1818037) q[2];
sx q[2];
rz(2.2369177) q[2];
rz(-0.71980643) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(-2.3864746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4633789) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(-0.71587193) q[0];
rz(2.8158358) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(-2.148927) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7870165) q[0];
sx q[0];
rz(-2.0873318) q[0];
sx q[0];
rz(1.2660962) q[0];
x q[1];
rz(2.2797212) q[2];
sx q[2];
rz(-2.3700691) q[2];
sx q[2];
rz(0.18415235) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0992972) q[1];
sx q[1];
rz(-2.9534833) q[1];
sx q[1];
rz(2.6204965) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1487001) q[3];
sx q[3];
rz(-1.143647) q[3];
sx q[3];
rz(-3.0689193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0830393) q[2];
sx q[2];
rz(-0.71648592) q[2];
sx q[2];
rz(-1.7129718) q[2];
rz(0.84609091) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(-1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5457526) q[0];
sx q[0];
rz(-1.4056453) q[0];
sx q[0];
rz(-0.76675057) q[0];
rz(1.2402361) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(-0.3516745) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2857916) q[0];
sx q[0];
rz(-1.8579146) q[0];
sx q[0];
rz(-2.48824) q[0];
rz(2.2064477) q[2];
sx q[2];
rz(-1.1543373) q[2];
sx q[2];
rz(0.93174975) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5778351) q[1];
sx q[1];
rz(-1.9655242) q[1];
sx q[1];
rz(1.5549591) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0760173) q[3];
sx q[3];
rz(-1.8540283) q[3];
sx q[3];
rz(-3.1056044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.19796431) q[2];
sx q[2];
rz(-1.6683656) q[2];
sx q[2];
rz(-0.53331214) q[2];
rz(2.7126281) q[3];
sx q[3];
rz(-0.52754378) q[3];
sx q[3];
rz(-2.0146446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-3.0241942) q[0];
sx q[0];
rz(-2.0485931) q[0];
sx q[0];
rz(0.54164106) q[0];
rz(-2.533124) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(-2.0419962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8655411) q[0];
sx q[0];
rz(-1.9380894) q[0];
sx q[0];
rz(2.0588576) q[0];
x q[1];
rz(1.0397644) q[2];
sx q[2];
rz(-1.4073939) q[2];
sx q[2];
rz(-2.642717) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6096566) q[1];
sx q[1];
rz(-1.7488639) q[1];
sx q[1];
rz(0.68105662) q[1];
rz(-pi) q[2];
rz(-2.8058365) q[3];
sx q[3];
rz(-0.81695518) q[3];
sx q[3];
rz(-2.9591065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5114484) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(-0.28298322) q[2];
rz(-0.7061559) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(-2.506536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.054984897) q[0];
sx q[0];
rz(-0.98751706) q[0];
sx q[0];
rz(1.5299861) q[0];
rz(0.64487547) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(-2.9842916) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0387715) q[0];
sx q[0];
rz(-2.4024995) q[0];
sx q[0];
rz(-0.13939136) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21643164) q[2];
sx q[2];
rz(-2.2955403) q[2];
sx q[2];
rz(-0.12491465) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0485222) q[1];
sx q[1];
rz(-1.375636) q[1];
sx q[1];
rz(-0.94961571) q[1];
rz(-0.39237202) q[3];
sx q[3];
rz(-2.4761768) q[3];
sx q[3];
rz(0.80442807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1823696) q[2];
sx q[2];
rz(-0.59252858) q[2];
sx q[2];
rz(2.3106993) q[2];
rz(-0.043878555) q[3];
sx q[3];
rz(-1.9079804) q[3];
sx q[3];
rz(-2.9390745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4956932) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(-0.010852531) q[0];
rz(-2.8649578) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(0.29327926) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0115259) q[0];
sx q[0];
rz(-1.5767326) q[0];
sx q[0];
rz(-0.54242155) q[0];
rz(-pi) q[1];
rz(2.9759334) q[2];
sx q[2];
rz(-2.207901) q[2];
sx q[2];
rz(0.13656244) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9353232) q[1];
sx q[1];
rz(-0.85201293) q[1];
sx q[1];
rz(1.000559) q[1];
rz(-pi) q[2];
rz(2.6777612) q[3];
sx q[3];
rz(-1.2004735) q[3];
sx q[3];
rz(0.90255373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.04348065) q[2];
sx q[2];
rz(-0.094823368) q[2];
sx q[2];
rz(0.088767178) q[2];
rz(2.8914715) q[3];
sx q[3];
rz(-2.0177896) q[3];
sx q[3];
rz(2.5626101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.0639169) q[0];
rz(0.44257277) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(-1.7907422) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8352855) q[0];
sx q[0];
rz(-1.4915691) q[0];
sx q[0];
rz(1.6552734) q[0];
rz(-pi) q[1];
rz(0.66019085) q[2];
sx q[2];
rz(-1.5056416) q[2];
sx q[2];
rz(2.8934663) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6752154) q[1];
sx q[1];
rz(-1.5870759) q[1];
sx q[1];
rz(2.0051763) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10339046) q[3];
sx q[3];
rz(-1.3500462) q[3];
sx q[3];
rz(3.0605928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0315447) q[2];
sx q[2];
rz(-1.2908547) q[2];
sx q[2];
rz(-2.6943977) q[2];
rz(-1.3859008) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0376301) q[0];
sx q[0];
rz(-2.5314405) q[0];
sx q[0];
rz(-2.2274341) q[0];
rz(0.42189235) q[2];
sx q[2];
rz(-1.5483529) q[2];
sx q[2];
rz(2.9107712) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1093724) q[1];
sx q[1];
rz(-1.6175555) q[1];
sx q[1];
rz(-2.1857775) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8923558) q[3];
sx q[3];
rz(-1.7064629) q[3];
sx q[3];
rz(0.78007573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.85049373) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(-2.4492241) q[2];
rz(-2.678357) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(-2.0326116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2904084) q[0];
sx q[0];
rz(-1.6155227) q[0];
sx q[0];
rz(0.68646705) q[0];
rz(2.6295173) q[1];
sx q[1];
rz(-2.8072186) q[1];
sx q[1];
rz(1.0288815) q[1];
rz(-0.26906536) q[2];
sx q[2];
rz(-1.27956) q[2];
sx q[2];
rz(0.051824311) q[2];
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
