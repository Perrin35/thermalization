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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.046145) q[0];
sx q[0];
rz(-2.5382222) q[0];
sx q[0];
rz(1.7741816) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6167294) q[2];
sx q[2];
rz(-2.7537991) q[2];
sx q[2];
rz(2.2847069) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9909518) q[1];
sx q[1];
rz(-0.71326637) q[1];
sx q[1];
rz(0.13891061) q[1];
rz(-2.9133137) q[3];
sx q[3];
rz(-2.7242959) q[3];
sx q[3];
rz(1.8018064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6266142) q[2];
sx q[2];
rz(-0.86003059) q[2];
sx q[2];
rz(2.853945) q[2];
rz(1.3927762) q[3];
sx q[3];
rz(-2.1578622) q[3];
sx q[3];
rz(-2.0387409) q[3];
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
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87067938) q[0];
sx q[0];
rz(-0.7551071) q[0];
sx q[0];
rz(2.5640008) q[0];
rz(1.5354935) q[1];
sx q[1];
rz(-2.0237193) q[1];
sx q[1];
rz(-1.5637406) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751511) q[0];
sx q[0];
rz(-1.6834714) q[0];
sx q[0];
rz(2.9215647) q[0];
rz(-pi) q[1];
rz(-1.2969442) q[2];
sx q[2];
rz(-1.8631862) q[2];
sx q[2];
rz(1.3016303) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4430122) q[1];
sx q[1];
rz(-0.46291446) q[1];
sx q[1];
rz(1.3742616) q[1];
rz(-0.79264499) q[3];
sx q[3];
rz(-2.9229197) q[3];
sx q[3];
rz(-1.6817026) q[3];
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
rz(-1.1632495) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(-0.64003402) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0796233) q[0];
sx q[0];
rz(-0.29456961) q[0];
sx q[0];
rz(0.18297718) q[0];
rz(0.043116365) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(-1.9972237) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9632918) q[0];
sx q[0];
rz(-0.68094567) q[0];
sx q[0];
rz(2.2292024) q[0];
x q[1];
rz(-1.3470596) q[2];
sx q[2];
rz(-1.9800595) q[2];
sx q[2];
rz(-0.70659107) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5444762) q[1];
sx q[1];
rz(-2.5529478) q[1];
sx q[1];
rz(2.4878923) q[1];
rz(-pi) q[2];
rz(-0.28762443) q[3];
sx q[3];
rz(-1.5737185) q[3];
sx q[3];
rz(2.6401816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.019505067) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(2.2369177) q[2];
rz(2.4217862) q[3];
sx q[3];
rz(-2.7147839) q[3];
sx q[3];
rz(2.3864746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.4633789) q[0];
sx q[0];
rz(-0.97020522) q[0];
sx q[0];
rz(2.4257207) q[0];
rz(-2.8158358) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(-2.148927) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0794605) q[0];
sx q[0];
rz(-1.3068763) q[0];
sx q[0];
rz(2.604565) q[0];
rz(0.86187141) q[2];
sx q[2];
rz(-2.3700691) q[2];
sx q[2];
rz(2.9574403) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6550733) q[1];
sx q[1];
rz(-1.7336978) q[1];
sx q[1];
rz(-1.4763114) q[1];
x q[2];
rz(-2.2654815) q[3];
sx q[3];
rz(-2.4377341) q[3];
sx q[3];
rz(-1.0775523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0830393) q[2];
sx q[2];
rz(-0.71648592) q[2];
sx q[2];
rz(1.4286208) q[2];
rz(-0.84609091) q[3];
sx q[3];
rz(-1.5139791) q[3];
sx q[3];
rz(1.9184453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5457526) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(2.3748421) q[0];
rz(1.9013566) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(0.3516745) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2131166) q[0];
sx q[0];
rz(-0.94841829) q[0];
sx q[0];
rz(-1.9268131) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93514498) q[2];
sx q[2];
rz(-1.1543373) q[2];
sx q[2];
rz(0.93174975) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60492086) q[1];
sx q[1];
rz(-2.7465638) q[1];
sx q[1];
rz(-3.1035963) q[1];
x q[2];
rz(2.8204927) q[3];
sx q[3];
rz(-1.0874815) q[3];
sx q[3];
rz(-1.6881642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19796431) q[2];
sx q[2];
rz(-1.4732271) q[2];
sx q[2];
rz(0.53331214) q[2];
rz(2.7126281) q[3];
sx q[3];
rz(-2.6140489) q[3];
sx q[3];
rz(-1.126948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0241942) q[0];
sx q[0];
rz(-2.0485931) q[0];
sx q[0];
rz(2.5999516) q[0];
rz(-0.60846865) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(2.0419962) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8894419) q[0];
sx q[0];
rz(-0.60177207) q[0];
sx q[0];
rz(-2.2579231) q[0];
rz(-pi) q[1];
rz(1.8855368) q[2];
sx q[2];
rz(-2.5883) q[2];
sx q[2];
rz(1.7994583) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.82367831) q[1];
sx q[1];
rz(-2.4412529) q[1];
sx q[1];
rz(0.27842303) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3533456) q[3];
sx q[3];
rz(-1.3282093) q[3];
sx q[3];
rz(1.6227674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
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
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054984897) q[0];
sx q[0];
rz(-0.98751706) q[0];
sx q[0];
rz(1.5299861) q[0];
rz(-0.64487547) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(-0.15730102) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3646506) q[0];
sx q[0];
rz(-1.4770664) q[0];
sx q[0];
rz(-0.73424299) q[0];
x q[1];
rz(-2.925161) q[2];
sx q[2];
rz(-0.84605233) q[2];
sx q[2];
rz(-0.12491465) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3992608) q[1];
sx q[1];
rz(-0.6472339) q[1];
sx q[1];
rz(-1.8982235) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7492206) q[3];
sx q[3];
rz(-2.4761768) q[3];
sx q[3];
rz(-2.3371646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9592231) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(2.3106993) q[2];
rz(-0.043878555) q[3];
sx q[3];
rz(-1.9079804) q[3];
sx q[3];
rz(0.20251814) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4956932) q[0];
sx q[0];
rz(-0.93911397) q[0];
sx q[0];
rz(3.1307401) q[0];
rz(-2.8649578) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(0.29327926) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1300668) q[0];
sx q[0];
rz(-1.5648601) q[0];
sx q[0];
rz(0.54242155) q[0];
rz(-0.16565928) q[2];
sx q[2];
rz(-0.93369166) q[2];
sx q[2];
rz(3.0050302) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3775023) q[1];
sx q[1];
rz(-1.1524156) q[1];
sx q[1];
rz(2.3368895) q[1];
rz(-pi) q[2];
rz(-1.1612438) q[3];
sx q[3];
rz(-1.1405986) q[3];
sx q[3];
rz(-0.48914117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.04348065) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(0.088767178) q[2];
rz(2.8914715) q[3];
sx q[3];
rz(-1.1238031) q[3];
sx q[3];
rz(-2.5626101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1254743) q[0];
sx q[0];
rz(-3.0258412) q[0];
sx q[0];
rz(2.3257757) q[0];
x q[1];
rz(-0.10599372) q[2];
sx q[2];
rz(-0.66291891) q[2];
sx q[2];
rz(-1.4063327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.46637725) q[1];
sx q[1];
rz(-1.5870759) q[1];
sx q[1];
rz(2.0051763) q[1];
rz(-1.1397347) q[3];
sx q[3];
rz(-0.24340478) q[3];
sx q[3];
rz(0.36153015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0315447) q[2];
sx q[2];
rz(-1.2908547) q[2];
sx q[2];
rz(-2.6943977) q[2];
rz(1.7556919) q[3];
sx q[3];
rz(-0.30062506) q[3];
sx q[3];
rz(-0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31496012) q[0];
sx q[0];
rz(-1.5216014) q[0];
sx q[0];
rz(-0.78053027) q[0];
rz(-2.7138117) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(2.9719877) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34925941) q[0];
sx q[0];
rz(-2.0418641) q[0];
sx q[0];
rz(2.7381998) q[0];
rz(-pi) q[1];
rz(1.595396) q[2];
sx q[2];
rz(-1.9925756) q[2];
sx q[2];
rz(-1.8116902) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6360213) q[1];
sx q[1];
rz(-2.1850056) q[1];
sx q[1];
rz(-3.0843656) q[1];
x q[2];
rz(-2.9986936) q[3];
sx q[3];
rz(-1.8892965) q[3];
sx q[3];
rz(2.3058476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.85049373) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(-2.4492241) q[2];
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
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2904084) q[0];
sx q[0];
rz(-1.52607) q[0];
sx q[0];
rz(-2.4551256) q[0];
rz(0.51207536) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(-1.2693263) q[2];
sx q[2];
rz(-1.3133247) q[2];
sx q[2];
rz(1.7016344) q[2];
rz(1.3626171) q[3];
sx q[3];
rz(-1.4559742) q[3];
sx q[3];
rz(-2.5696587) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
