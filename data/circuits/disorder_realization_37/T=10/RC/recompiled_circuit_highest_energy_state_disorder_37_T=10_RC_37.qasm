OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3353424) q[0];
sx q[0];
rz(3.5080533) q[0];
sx q[0];
rz(8.9656497) q[0];
rz(2.4899809) q[1];
sx q[1];
rz(-1.4067283) q[1];
sx q[1];
rz(-0.23174098) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36759672) q[0];
sx q[0];
rz(-1.2277506) q[0];
sx q[0];
rz(1.6438707) q[0];
rz(-pi) q[1];
rz(-2.5951573) q[2];
sx q[2];
rz(-1.6197512) q[2];
sx q[2];
rz(-0.10412439) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.79855483) q[1];
sx q[1];
rz(-2.3904368) q[1];
sx q[1];
rz(0.57111528) q[1];
x q[2];
rz(-0.11756331) q[3];
sx q[3];
rz(-0.34053482) q[3];
sx q[3];
rz(-1.7528319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0509402) q[2];
sx q[2];
rz(-2.6068164) q[2];
sx q[2];
rz(1.862662) q[2];
rz(-2.2517962) q[3];
sx q[3];
rz(-1.5225478) q[3];
sx q[3];
rz(-2.8029627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58562529) q[0];
sx q[0];
rz(-2.2323759) q[0];
sx q[0];
rz(2.2157748) q[0];
rz(1.0649118) q[1];
sx q[1];
rz(-1.5644667) q[1];
sx q[1];
rz(0.085478641) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8877713) q[0];
sx q[0];
rz(-1.5740875) q[0];
sx q[0];
rz(-0.38981593) q[0];
rz(-pi) q[1];
rz(-0.95909848) q[2];
sx q[2];
rz(-1.8597721) q[2];
sx q[2];
rz(-2.0069569) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2457471) q[1];
sx q[1];
rz(-1.8008044) q[1];
sx q[1];
rz(1.0279015) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68686266) q[3];
sx q[3];
rz(-1.1986889) q[3];
sx q[3];
rz(2.5290979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.815879) q[2];
sx q[2];
rz(-1.4789944) q[2];
sx q[2];
rz(2.7549287) q[2];
rz(1.2169085) q[3];
sx q[3];
rz(-0.34438008) q[3];
sx q[3];
rz(-2.9215422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81095186) q[0];
sx q[0];
rz(-2.7702259) q[0];
sx q[0];
rz(-2.6445342) q[0];
rz(-1.0423543) q[1];
sx q[1];
rz(-0.94756871) q[1];
sx q[1];
rz(0.29464468) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9884856) q[0];
sx q[0];
rz(-1.238958) q[0];
sx q[0];
rz(1.6436623) q[0];
rz(-pi) q[1];
rz(-0.23068409) q[2];
sx q[2];
rz(-0.65197021) q[2];
sx q[2];
rz(-2.1459652) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4528447) q[1];
sx q[1];
rz(-2.0032681) q[1];
sx q[1];
rz(2.2562863) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1717242) q[3];
sx q[3];
rz(-1.4517541) q[3];
sx q[3];
rz(2.5210019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7872539) q[2];
sx q[2];
rz(-2.2806809) q[2];
sx q[2];
rz(-1.5798689) q[2];
rz(-2.0653557) q[3];
sx q[3];
rz(-1.302224) q[3];
sx q[3];
rz(2.2126183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8266066) q[0];
sx q[0];
rz(-1.5262693) q[0];
sx q[0];
rz(2.9122747) q[0];
rz(1.6353105) q[1];
sx q[1];
rz(-1.6835667) q[1];
sx q[1];
rz(0.062072676) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44767013) q[0];
sx q[0];
rz(-1.0167443) q[0];
sx q[0];
rz(-0.5918795) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96458413) q[2];
sx q[2];
rz(-0.46614376) q[2];
sx q[2];
rz(0.073748253) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5531299) q[1];
sx q[1];
rz(-2.3847918) q[1];
sx q[1];
rz(-2.7314145) q[1];
rz(-pi) q[2];
rz(-1.4775949) q[3];
sx q[3];
rz(-1.3706932) q[3];
sx q[3];
rz(-0.72678185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.092992358) q[2];
sx q[2];
rz(-2.2255662) q[2];
sx q[2];
rz(-1.3108866) q[2];
rz(3.1033031) q[3];
sx q[3];
rz(-1.7891276) q[3];
sx q[3];
rz(0.37461764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6467317) q[0];
sx q[0];
rz(-0.72135389) q[0];
sx q[0];
rz(0.89299655) q[0];
rz(-2.112174) q[1];
sx q[1];
rz(-2.4118377) q[1];
sx q[1];
rz(-3.0457048) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9150118) q[0];
sx q[0];
rz(-1.1553197) q[0];
sx q[0];
rz(1.0092495) q[0];
rz(-0.070438373) q[2];
sx q[2];
rz(-1.6595006) q[2];
sx q[2];
rz(0.34102893) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3600625) q[1];
sx q[1];
rz(-1.9340098) q[1];
sx q[1];
rz(-2.5468154) q[1];
x q[2];
rz(0.60378243) q[3];
sx q[3];
rz(-0.60294916) q[3];
sx q[3];
rz(-0.83151885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15478495) q[2];
sx q[2];
rz(-1.751535) q[2];
sx q[2];
rz(0.42547697) q[2];
rz(0.42243877) q[3];
sx q[3];
rz(-2.0368302) q[3];
sx q[3];
rz(-2.6384242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7285889) q[0];
sx q[0];
rz(-0.63218963) q[0];
sx q[0];
rz(3.0244306) q[0];
rz(1.6475742) q[1];
sx q[1];
rz(-0.90227503) q[1];
sx q[1];
rz(-2.07043) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8534096) q[0];
sx q[0];
rz(-0.57153915) q[0];
sx q[0];
rz(-2.6009212) q[0];
rz(-2.5994634) q[2];
sx q[2];
rz(-1.8111992) q[2];
sx q[2];
rz(-1.2803276) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6630604) q[1];
sx q[1];
rz(-0.89420616) q[1];
sx q[1];
rz(0.29037906) q[1];
x q[2];
rz(-0.31757327) q[3];
sx q[3];
rz(-1.4151207) q[3];
sx q[3];
rz(2.6359216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5693207) q[2];
sx q[2];
rz(-2.5666777) q[2];
sx q[2];
rz(-3.102109) q[2];
rz(3.0651921) q[3];
sx q[3];
rz(-1.1698086) q[3];
sx q[3];
rz(0.45342818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4385248) q[0];
sx q[0];
rz(-2.330133) q[0];
sx q[0];
rz(-0.95034289) q[0];
rz(-1.9901265) q[1];
sx q[1];
rz(-1.5299503) q[1];
sx q[1];
rz(1.8340402) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81539736) q[0];
sx q[0];
rz(-1.895825) q[0];
sx q[0];
rz(0.23120489) q[0];
rz(-pi) q[1];
rz(-1.4979532) q[2];
sx q[2];
rz(-1.1955639) q[2];
sx q[2];
rz(0.90922395) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.79297355) q[1];
sx q[1];
rz(-1.359531) q[1];
sx q[1];
rz(-1.068685) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3780955) q[3];
sx q[3];
rz(-0.61470882) q[3];
sx q[3];
rz(2.4845882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.82039708) q[2];
sx q[2];
rz(-2.4615007) q[2];
sx q[2];
rz(-0.55366984) q[2];
rz(-0.6959483) q[3];
sx q[3];
rz(-1.4724933) q[3];
sx q[3];
rz(-1.6863916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11947908) q[0];
sx q[0];
rz(-1.6895634) q[0];
sx q[0];
rz(-0.30690646) q[0];
rz(-1.2762997) q[1];
sx q[1];
rz(-0.46638322) q[1];
sx q[1];
rz(1.2409522) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0690445) q[0];
sx q[0];
rz(-2.6676151) q[0];
sx q[0];
rz(1.9182253) q[0];
x q[1];
rz(2.5503301) q[2];
sx q[2];
rz(-2.4545672) q[2];
sx q[2];
rz(0.20461539) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7684264) q[1];
sx q[1];
rz(-1.8796928) q[1];
sx q[1];
rz(-1.7684446) q[1];
rz(-pi) q[2];
rz(-0.54540619) q[3];
sx q[3];
rz(-2.4562803) q[3];
sx q[3];
rz(-2.3092676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.79157311) q[2];
sx q[2];
rz(-2.3363523) q[2];
sx q[2];
rz(1.8512858) q[2];
rz(0.6811412) q[3];
sx q[3];
rz(-2.2414424) q[3];
sx q[3];
rz(1.6397938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8660368) q[0];
sx q[0];
rz(-2.10502) q[0];
sx q[0];
rz(1.9679029) q[0];
rz(-2.5545919) q[1];
sx q[1];
rz(-2.360011) q[1];
sx q[1];
rz(-0.079708286) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22686887) q[0];
sx q[0];
rz(-2.7023315) q[0];
sx q[0];
rz(-2.125581) q[0];
x q[1];
rz(2.6775421) q[2];
sx q[2];
rz(-2.8949605) q[2];
sx q[2];
rz(0.76619785) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.6882872) q[1];
sx q[1];
rz(-0.88677471) q[1];
sx q[1];
rz(0.75835336) q[1];
rz(-0.80612889) q[3];
sx q[3];
rz(-1.3622096) q[3];
sx q[3];
rz(-0.94193469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54032636) q[2];
sx q[2];
rz(-1.2286295) q[2];
sx q[2];
rz(-1.3339174) q[2];
rz(-1.5506844) q[3];
sx q[3];
rz(-2.1315137) q[3];
sx q[3];
rz(0.18868119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.659336) q[0];
sx q[0];
rz(-1.5631258) q[0];
sx q[0];
rz(-1.851086) q[0];
rz(3.105063) q[1];
sx q[1];
rz(-1.9772269) q[1];
sx q[1];
rz(-2.0972924) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.288702) q[0];
sx q[0];
rz(-1.9963309) q[0];
sx q[0];
rz(2.8136926) q[0];
rz(2.2453111) q[2];
sx q[2];
rz(-2.7283629) q[2];
sx q[2];
rz(-0.21597029) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.579291) q[1];
sx q[1];
rz(-2.3962483) q[1];
sx q[1];
rz(-1.7991123) q[1];
rz(2.815991) q[3];
sx q[3];
rz(-0.34593098) q[3];
sx q[3];
rz(2.9671362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2269939) q[2];
sx q[2];
rz(-1.031216) q[2];
sx q[2];
rz(-1.3247789) q[2];
rz(-2.3850208) q[3];
sx q[3];
rz(-2.2533267) q[3];
sx q[3];
rz(-0.85048401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4751547) q[0];
sx q[0];
rz(-1.7378687) q[0];
sx q[0];
rz(0.73874656) q[0];
rz(-1.7186164) q[1];
sx q[1];
rz(-1.1971133) q[1];
sx q[1];
rz(-2.8308629) q[1];
rz(-2.4074211) q[2];
sx q[2];
rz(-1.894507) q[2];
sx q[2];
rz(1.4121216) q[2];
rz(1.1005836) q[3];
sx q[3];
rz(-0.96405021) q[3];
sx q[3];
rz(-1.8586803) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
