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
rz(0.75594354) q[0];
sx q[0];
rz(-2.8109901) q[0];
sx q[0];
rz(-0.27275738) q[0];
rz(2.7859712) q[1];
sx q[1];
rz(-2.1115117) q[1];
sx q[1];
rz(-1.570809) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8518675) q[0];
sx q[0];
rz(-3.0209732) q[0];
sx q[0];
rz(-2.2096746) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5982389) q[2];
sx q[2];
rz(-2.5961868) q[2];
sx q[2];
rz(2.0871833) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.194549) q[1];
sx q[1];
rz(-1.5284555) q[1];
sx q[1];
rz(-0.6794828) q[1];
rz(-pi) q[2];
rz(0.36873534) q[3];
sx q[3];
rz(-1.6622541) q[3];
sx q[3];
rz(-1.9029332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.97751272) q[2];
sx q[2];
rz(-1.432632) q[2];
sx q[2];
rz(2.254503) q[2];
rz(-1.8767493) q[3];
sx q[3];
rz(-1.0586459) q[3];
sx q[3];
rz(0.023887159) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3905268) q[0];
sx q[0];
rz(-0.85705119) q[0];
sx q[0];
rz(0.29399011) q[0];
rz(1.7901621) q[1];
sx q[1];
rz(-2.0452812) q[1];
sx q[1];
rz(-1.9757804) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6751926) q[0];
sx q[0];
rz(-1.2443911) q[0];
sx q[0];
rz(-0.2194039) q[0];
rz(-pi) q[1];
rz(1.9306669) q[2];
sx q[2];
rz(-3.0635298) q[2];
sx q[2];
rz(0.087457267) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.15951751) q[1];
sx q[1];
rz(-1.0613835) q[1];
sx q[1];
rz(-2.0768836) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8610624) q[3];
sx q[3];
rz(-2.120388) q[3];
sx q[3];
rz(-2.275718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4531648) q[2];
sx q[2];
rz(-2.5025949) q[2];
sx q[2];
rz(1.2028018) q[2];
rz(1.5419675) q[3];
sx q[3];
rz(-1.3381693) q[3];
sx q[3];
rz(-0.28551027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-0.11005345) q[0];
sx q[0];
rz(-0.078131214) q[0];
sx q[0];
rz(-1.9670638) q[0];
rz(-2.166676) q[1];
sx q[1];
rz(-2.2623623) q[1];
sx q[1];
rz(2.6851795) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83041265) q[0];
sx q[0];
rz(-1.5390656) q[0];
sx q[0];
rz(1.8422442) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9554702) q[2];
sx q[2];
rz(-1.3403424) q[2];
sx q[2];
rz(1.3467033) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3756936) q[1];
sx q[1];
rz(-2.1296892) q[1];
sx q[1];
rz(0.8984579) q[1];
rz(-1.0393082) q[3];
sx q[3];
rz(-0.78954299) q[3];
sx q[3];
rz(-1.9122151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6156893) q[2];
sx q[2];
rz(-0.38280767) q[2];
sx q[2];
rz(2.3846386) q[2];
rz(3.1225539) q[3];
sx q[3];
rz(-1.8502539) q[3];
sx q[3];
rz(2.3058057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2794936) q[0];
sx q[0];
rz(-2.4308496) q[0];
sx q[0];
rz(-0.94974649) q[0];
rz(0.21385916) q[1];
sx q[1];
rz(-2.2016134) q[1];
sx q[1];
rz(-0.59923879) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1480162) q[0];
sx q[0];
rz(-2.2207119) q[0];
sx q[0];
rz(0.17493576) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9733359) q[2];
sx q[2];
rz(-2.0806542) q[2];
sx q[2];
rz(-2.5349701) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.487452) q[1];
sx q[1];
rz(-2.096744) q[1];
sx q[1];
rz(-0.24152813) q[1];
rz(-pi) q[2];
rz(-1.4324801) q[3];
sx q[3];
rz(-2.1475422) q[3];
sx q[3];
rz(-0.53046295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.63138258) q[2];
sx q[2];
rz(-1.6333132) q[2];
sx q[2];
rz(-2.317826) q[2];
rz(1.0907762) q[3];
sx q[3];
rz(-2.3407276) q[3];
sx q[3];
rz(-2.4763079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3983696) q[0];
sx q[0];
rz(-0.38341612) q[0];
sx q[0];
rz(-0.98633352) q[0];
rz(2.8979454) q[1];
sx q[1];
rz(-1.2657093) q[1];
sx q[1];
rz(-2.9820014) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8411113) q[0];
sx q[0];
rz(-2.2694025) q[0];
sx q[0];
rz(2.7270856) q[0];
rz(-0.56964376) q[2];
sx q[2];
rz(-2.2692371) q[2];
sx q[2];
rz(1.4739571) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1732542) q[1];
sx q[1];
rz(-0.45351492) q[1];
sx q[1];
rz(-1.3682925) q[1];
x q[2];
rz(-0.87017228) q[3];
sx q[3];
rz(-1.7990094) q[3];
sx q[3];
rz(0.96503497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39449447) q[2];
sx q[2];
rz(-1.2622086) q[2];
sx q[2];
rz(2.0060284) q[2];
rz(-0.46843946) q[3];
sx q[3];
rz(-1.2218385) q[3];
sx q[3];
rz(-0.92456094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.121345) q[0];
sx q[0];
rz(-2.0609042) q[0];
sx q[0];
rz(-2.8998846) q[0];
rz(-1.3555591) q[1];
sx q[1];
rz(-0.87239289) q[1];
sx q[1];
rz(2.2579069) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6304105) q[0];
sx q[0];
rz(-2.1616461) q[0];
sx q[0];
rz(0.972959) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68285791) q[2];
sx q[2];
rz(-2.747263) q[2];
sx q[2];
rz(-0.71428052) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.984772) q[1];
sx q[1];
rz(-1.3756911) q[1];
sx q[1];
rz(2.0181993) q[1];
rz(0.00022907654) q[3];
sx q[3];
rz(-1.159824) q[3];
sx q[3];
rz(-0.64490333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9200865) q[2];
sx q[2];
rz(-2.2111427) q[2];
sx q[2];
rz(2.7152756) q[2];
rz(2.2250371) q[3];
sx q[3];
rz(-1.5896268) q[3];
sx q[3];
rz(-2.0385273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78974462) q[0];
sx q[0];
rz(-1.177657) q[0];
sx q[0];
rz(1.1258997) q[0];
rz(-0.060983505) q[1];
sx q[1];
rz(-0.95044249) q[1];
sx q[1];
rz(1.0152063) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7506152) q[0];
sx q[0];
rz(-0.96834194) q[0];
sx q[0];
rz(1.0573122) q[0];
rz(-1.9092375) q[2];
sx q[2];
rz(-1.5012174) q[2];
sx q[2];
rz(1.4078946) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.2992933) q[1];
sx q[1];
rz(-1.1720578) q[1];
sx q[1];
rz(-0.066940424) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0281248) q[3];
sx q[3];
rz(-1.610329) q[3];
sx q[3];
rz(-3.0276445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8288377) q[2];
sx q[2];
rz(-1.2042896) q[2];
sx q[2];
rz(-2.915536) q[2];
rz(-2.1894646) q[3];
sx q[3];
rz(-3.003037) q[3];
sx q[3];
rz(2.3329195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1299745) q[0];
sx q[0];
rz(-1.3148774) q[0];
sx q[0];
rz(-0.34039482) q[0];
rz(0.099523425) q[1];
sx q[1];
rz(-2.0390022) q[1];
sx q[1];
rz(-1.5955101) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81170245) q[0];
sx q[0];
rz(-1.8970338) q[0];
sx q[0];
rz(-2.5094112) q[0];
rz(-1.7695707) q[2];
sx q[2];
rz(-0.27789657) q[2];
sx q[2];
rz(-0.0065491876) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9118285) q[1];
sx q[1];
rz(-0.10783261) q[1];
sx q[1];
rz(2.7417438) q[1];
rz(-pi) q[2];
rz(0.45507064) q[3];
sx q[3];
rz(-1.651147) q[3];
sx q[3];
rz(-0.72380607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16326363) q[2];
sx q[2];
rz(-1.949387) q[2];
sx q[2];
rz(2.1400129) q[2];
rz(-1.7892276) q[3];
sx q[3];
rz(-2.6632301) q[3];
sx q[3];
rz(-1.9549687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.96420646) q[0];
sx q[0];
rz(-2.0851676) q[0];
sx q[0];
rz(-0.00038432234) q[0];
rz(-2.7035825) q[1];
sx q[1];
rz(-1.6338467) q[1];
sx q[1];
rz(2.6843574) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9033636) q[0];
sx q[0];
rz(-1.6203124) q[0];
sx q[0];
rz(1.7286506) q[0];
rz(-pi) q[1];
rz(0.57827823) q[2];
sx q[2];
rz(-2.1077029) q[2];
sx q[2];
rz(1.121305) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.776062) q[1];
sx q[1];
rz(-1.9161822) q[1];
sx q[1];
rz(-1.3341435) q[1];
x q[2];
rz(-1.2401388) q[3];
sx q[3];
rz(-2.0539502) q[3];
sx q[3];
rz(-1.5548979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.75927258) q[2];
sx q[2];
rz(-2.1542532) q[2];
sx q[2];
rz(1.5116723) q[2];
rz(-0.090653732) q[3];
sx q[3];
rz(-1.0414618) q[3];
sx q[3];
rz(0.67971984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35336211) q[0];
sx q[0];
rz(-1.7857977) q[0];
sx q[0];
rz(2.8613388) q[0];
rz(1.7578112) q[1];
sx q[1];
rz(-0.46418142) q[1];
sx q[1];
rz(1.9162477) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89577755) q[0];
sx q[0];
rz(-1.4207463) q[0];
sx q[0];
rz(1.4729985) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2544075) q[2];
sx q[2];
rz(-1.2204183) q[2];
sx q[2];
rz(-0.35752192) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.16364747) q[1];
sx q[1];
rz(-2.3446348) q[1];
sx q[1];
rz(-2.7989788) q[1];
rz(-pi) q[2];
rz(2.8392196) q[3];
sx q[3];
rz(-1.8480715) q[3];
sx q[3];
rz(-2.759861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19259351) q[2];
sx q[2];
rz(-1.6152363) q[2];
sx q[2];
rz(2.4601649) q[2];
rz(-1.1000018) q[3];
sx q[3];
rz(-0.93183485) q[3];
sx q[3];
rz(-1.3070235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6901907) q[0];
sx q[0];
rz(-2.207353) q[0];
sx q[0];
rz(1.9718715) q[0];
rz(-2.7611217) q[1];
sx q[1];
rz(-2.4741551) q[1];
sx q[1];
rz(-2.9173775) q[1];
rz(0.26421122) q[2];
sx q[2];
rz(-0.88709863) q[2];
sx q[2];
rz(-1.8567793) q[2];
rz(2.1721645) q[3];
sx q[3];
rz(-0.68404861) q[3];
sx q[3];
rz(0.28290471) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
