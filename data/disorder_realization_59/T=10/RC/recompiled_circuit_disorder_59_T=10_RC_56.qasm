OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(2.3595915) q[0];
sx q[0];
rz(10.695988) q[0];
rz(0.27702364) q[1];
sx q[1];
rz(-0.47200534) q[1];
sx q[1];
rz(0.0013874887) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.815925) q[0];
sx q[0];
rz(-2.7357833) q[0];
sx q[0];
rz(-1.9947467) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0543047) q[2];
sx q[2];
rz(-2.6929571) q[2];
sx q[2];
rz(1.0686312) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.672294) q[1];
sx q[1];
rz(-2.0672332) q[1];
sx q[1];
rz(-1.6507571) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1625159) q[3];
sx q[3];
rz(-0.43833971) q[3];
sx q[3];
rz(-1.0867659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9871621) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(0.74938613) q[2];
rz(-2.1253712) q[3];
sx q[3];
rz(-1.1775492) q[3];
sx q[3];
rz(2.7367676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4352903) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(-2.170927) q[0];
rz(-2.1043815) q[1];
sx q[1];
rz(-1.7036006) q[1];
sx q[1];
rz(2.326139) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97074189) q[0];
sx q[0];
rz(-2.4903957) q[0];
sx q[0];
rz(-0.09911508) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.843812) q[2];
sx q[2];
rz(-0.85999876) q[2];
sx q[2];
rz(-0.057991512) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7533469) q[1];
sx q[1];
rz(-1.2699632) q[1];
sx q[1];
rz(-1.6080329) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8498597) q[3];
sx q[3];
rz(-2.4749304) q[3];
sx q[3];
rz(0.093689703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4619535) q[2];
sx q[2];
rz(-1.5660428) q[2];
sx q[2];
rz(-0.63278502) q[2];
rz(1.1535545) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(-0.30953428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84045029) q[0];
sx q[0];
rz(-1.2121032) q[0];
sx q[0];
rz(-2.2667623) q[0];
rz(1.3300928) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(0.99951807) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32854983) q[0];
sx q[0];
rz(-1.8694436) q[0];
sx q[0];
rz(0.15655984) q[0];
rz(-pi) q[1];
rz(2.2314084) q[2];
sx q[2];
rz(-1.0599469) q[2];
sx q[2];
rz(0.78188932) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1332902) q[1];
sx q[1];
rz(-1.8255207) q[1];
sx q[1];
rz(-0.50526527) q[1];
x q[2];
rz(-2.0796892) q[3];
sx q[3];
rz(-1.6403927) q[3];
sx q[3];
rz(0.059046179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6040566) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(0.58004722) q[2];
rz(0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(-1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.375305) q[0];
sx q[0];
rz(-1.5505318) q[0];
sx q[0];
rz(0.91039175) q[0];
rz(-0.45122775) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(-2.8667563) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66514689) q[0];
sx q[0];
rz(-3.0229212) q[0];
sx q[0];
rz(0.84400405) q[0];
rz(-pi) q[1];
rz(-1.1103815) q[2];
sx q[2];
rz(-1.5722256) q[2];
sx q[2];
rz(1.9542076) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.30157629) q[1];
sx q[1];
rz(-2.3641572) q[1];
sx q[1];
rz(2.2393054) q[1];
x q[2];
rz(0.87423012) q[3];
sx q[3];
rz(-1.6682373) q[3];
sx q[3];
rz(1.8271354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.794902) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(1.5083195) q[2];
rz(-1.1446965) q[3];
sx q[3];
rz(-0.73994023) q[3];
sx q[3];
rz(-2.9798853) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65790025) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(-0.99779469) q[0];
rz(-2.9580341) q[1];
sx q[1];
rz(-1.654637) q[1];
sx q[1];
rz(1.6246187) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8109587) q[0];
sx q[0];
rz(-1.9348382) q[0];
sx q[0];
rz(0.90322687) q[0];
rz(-pi) q[1];
rz(-0.60913182) q[2];
sx q[2];
rz(-0.62437781) q[2];
sx q[2];
rz(2.6513211) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6301873) q[1];
sx q[1];
rz(-1.9747707) q[1];
sx q[1];
rz(-1.1955111) q[1];
rz(-0.59668221) q[3];
sx q[3];
rz(-1.0718857) q[3];
sx q[3];
rz(-0.98017207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3395485) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(2.8175763) q[2];
rz(-1.3230532) q[3];
sx q[3];
rz(-0.75606212) q[3];
sx q[3];
rz(-1.5312622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762387) q[0];
sx q[0];
rz(-2.1026251) q[0];
sx q[0];
rz(-1.2639686) q[0];
rz(0.91066796) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(-2.8009159) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7621988) q[0];
sx q[0];
rz(-1.5389991) q[0];
sx q[0];
rz(-1.5048774) q[0];
rz(2.3855626) q[2];
sx q[2];
rz(-1.7871734) q[2];
sx q[2];
rz(1.1083958) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1813982) q[1];
sx q[1];
rz(-1.3409233) q[1];
sx q[1];
rz(-0.48744907) q[1];
rz(2.7192781) q[3];
sx q[3];
rz(-1.1186244) q[3];
sx q[3];
rz(-0.57002588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.95057758) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(0.77159709) q[2];
rz(-0.54780444) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(2.5781393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9724378) q[0];
sx q[0];
rz(-1.6945524) q[0];
sx q[0];
rz(-0.34564885) q[0];
rz(3.0787643) q[1];
sx q[1];
rz(-0.47880104) q[1];
sx q[1];
rz(2.6766434) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9879887) q[0];
sx q[0];
rz(-2.0360332) q[0];
sx q[0];
rz(0.21501712) q[0];
rz(-0.49352383) q[2];
sx q[2];
rz(-2.0888622) q[2];
sx q[2];
rz(-1.354419) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.17864922) q[1];
sx q[1];
rz(-1.6766251) q[1];
sx q[1];
rz(-0.21957285) q[1];
rz(-pi) q[2];
rz(-1.5246478) q[3];
sx q[3];
rz(-1.5966291) q[3];
sx q[3];
rz(-1.8125364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0039625) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(-2.8239992) q[2];
rz(-2.5701304) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6222318) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(-0.28433329) q[0];
rz(-0.55150664) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(3.0632339) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4975472) q[0];
sx q[0];
rz(-2.3345778) q[0];
sx q[0];
rz(3.0456196) q[0];
x q[1];
rz(-0.80348357) q[2];
sx q[2];
rz(-0.33499559) q[2];
sx q[2];
rz(1.4575046) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.52787493) q[1];
sx q[1];
rz(-0.67786874) q[1];
sx q[1];
rz(1.8498969) q[1];
x q[2];
rz(1.7275229) q[3];
sx q[3];
rz(-1.6046451) q[3];
sx q[3];
rz(-1.0364929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7408961) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(-0.25137869) q[2];
rz(-2.5583983) q[3];
sx q[3];
rz(-1.1116894) q[3];
sx q[3];
rz(-1.4096227) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3437929) q[0];
sx q[0];
rz(-3.0594337) q[0];
sx q[0];
rz(0.051368512) q[0];
rz(-2.2180166) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(0.87402469) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4955935) q[0];
sx q[0];
rz(-2.3388303) q[0];
sx q[0];
rz(1.5828703) q[0];
rz(-pi) q[1];
rz(2.2888695) q[2];
sx q[2];
rz(-0.52041473) q[2];
sx q[2];
rz(-0.66327099) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.28305001) q[1];
sx q[1];
rz(-2.0701323) q[1];
sx q[1];
rz(-3.0602171) q[1];
rz(-pi) q[2];
rz(0.47253982) q[3];
sx q[3];
rz(-0.66415411) q[3];
sx q[3];
rz(3.1411375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41436568) q[2];
sx q[2];
rz(-2.3870769) q[2];
sx q[2];
rz(0.61974636) q[2];
rz(1.9571346) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(1.363389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1353564) q[0];
sx q[0];
rz(-1.0422491) q[0];
sx q[0];
rz(2.4172879) q[0];
rz(-2.9528217) q[1];
sx q[1];
rz(-0.17938463) q[1];
sx q[1];
rz(1.1788517) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41335426) q[0];
sx q[0];
rz(-1.2984707) q[0];
sx q[0];
rz(-2.9288835) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63209052) q[2];
sx q[2];
rz(-0.1569911) q[2];
sx q[2];
rz(2.0096411) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9465543) q[1];
sx q[1];
rz(-2.185501) q[1];
sx q[1];
rz(-0.22919319) q[1];
x q[2];
rz(1.7970656) q[3];
sx q[3];
rz(-1.0591649) q[3];
sx q[3];
rz(-0.98480485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0835691) q[2];
sx q[2];
rz(-2.100779) q[2];
sx q[2];
rz(-0.89938346) q[2];
rz(0.87456885) q[3];
sx q[3];
rz(-2.7159297) q[3];
sx q[3];
rz(-1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.62109229) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(2.3836366) q[1];
sx q[1];
rz(-1.4823722) q[1];
sx q[1];
rz(-1.4588251) q[1];
rz(-1.6451251) q[2];
sx q[2];
rz(-0.44114124) q[2];
sx q[2];
rz(0.91976358) q[2];
rz(-0.48537985) q[3];
sx q[3];
rz(-2.831922) q[3];
sx q[3];
rz(0.61556863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
