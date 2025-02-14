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
rz(3.0313015) q[0];
sx q[0];
rz(-1.8561441) q[0];
sx q[0];
rz(0.31874803) q[0];
rz(1.1822074) q[1];
sx q[1];
rz(-1.6637586) q[1];
sx q[1];
rz(-1.1629265) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4773552) q[0];
sx q[0];
rz(-1.2384982) q[0];
sx q[0];
rz(-2.8927781) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0672795) q[2];
sx q[2];
rz(-0.19467672) q[2];
sx q[2];
rz(-3.1271324) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8713177) q[1];
sx q[1];
rz(-1.8757038) q[1];
sx q[1];
rz(-0.31590806) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5507023) q[3];
sx q[3];
rz(-2.1455975) q[3];
sx q[3];
rz(1.8679096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9691539) q[2];
sx q[2];
rz(-2.6814851) q[2];
sx q[2];
rz(-0.13574204) q[2];
rz(2.8067449) q[3];
sx q[3];
rz(-1.9183153) q[3];
sx q[3];
rz(2.7443583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927476) q[0];
sx q[0];
rz(-2.1362342) q[0];
sx q[0];
rz(0.94648615) q[0];
rz(-1.1874366) q[1];
sx q[1];
rz(-2.5577736) q[1];
sx q[1];
rz(-2.8576287) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89881247) q[0];
sx q[0];
rz(-1.3889342) q[0];
sx q[0];
rz(-2.9480431) q[0];
rz(-pi) q[1];
rz(-2.9049529) q[2];
sx q[2];
rz(-2.2556117) q[2];
sx q[2];
rz(-0.94517148) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3888106) q[1];
sx q[1];
rz(-1.4335263) q[1];
sx q[1];
rz(-2.7146882) q[1];
rz(-pi) q[2];
rz(-0.3881298) q[3];
sx q[3];
rz(-1.8947269) q[3];
sx q[3];
rz(-2.7517954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2076608) q[2];
sx q[2];
rz(-1.2958823) q[2];
sx q[2];
rz(2.3409823) q[2];
rz(-1.8396395) q[3];
sx q[3];
rz(-2.0079565) q[3];
sx q[3];
rz(-3.083526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0170853) q[0];
sx q[0];
rz(-0.43211102) q[0];
sx q[0];
rz(2.3486163) q[0];
rz(2.4555581) q[1];
sx q[1];
rz(-1.425309) q[1];
sx q[1];
rz(2.3763903) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6128639) q[0];
sx q[0];
rz(-2.9224797) q[0];
sx q[0];
rz(1.8414453) q[0];
rz(0.051570895) q[2];
sx q[2];
rz(-1.5386536) q[2];
sx q[2];
rz(0.78173897) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6531892) q[1];
sx q[1];
rz(-0.81522664) q[1];
sx q[1];
rz(1.9550187) q[1];
rz(-pi) q[2];
rz(0.92242494) q[3];
sx q[3];
rz(-2.2545345) q[3];
sx q[3];
rz(0.91372638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5746295) q[2];
sx q[2];
rz(-2.4506863) q[2];
sx q[2];
rz(1.2464657) q[2];
rz(1.9555107) q[3];
sx q[3];
rz(-1.3776774) q[3];
sx q[3];
rz(-0.20868364) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74801159) q[0];
sx q[0];
rz(-1.8210541) q[0];
sx q[0];
rz(3.0082974) q[0];
rz(2.8721299) q[1];
sx q[1];
rz(-0.91454426) q[1];
sx q[1];
rz(-2.6752313) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5488777) q[0];
sx q[0];
rz(-0.60564954) q[0];
sx q[0];
rz(2.9306508) q[0];
x q[1];
rz(0.94302098) q[2];
sx q[2];
rz(-1.4157214) q[2];
sx q[2];
rz(2.6076406) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7940841) q[1];
sx q[1];
rz(-1.7219436) q[1];
sx q[1];
rz(-2.529649) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6367988) q[3];
sx q[3];
rz(-1.0579234) q[3];
sx q[3];
rz(-0.51067019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.93126297) q[2];
sx q[2];
rz(-0.42079058) q[2];
sx q[2];
rz(-0.24410625) q[2];
rz(-1.3611475) q[3];
sx q[3];
rz(-0.94435349) q[3];
sx q[3];
rz(-1.2066427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4162783) q[0];
sx q[0];
rz(-0.84753528) q[0];
sx q[0];
rz(-0.086294802) q[0];
rz(1.3823973) q[1];
sx q[1];
rz(-2.081213) q[1];
sx q[1];
rz(-1.8550526) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70885284) q[0];
sx q[0];
rz(-1.7403495) q[0];
sx q[0];
rz(1.1772593) q[0];
rz(2.3628144) q[2];
sx q[2];
rz(-1.4113103) q[2];
sx q[2];
rz(-0.82568141) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.55505048) q[1];
sx q[1];
rz(-2.4355531) q[1];
sx q[1];
rz(-1.4439739) q[1];
x q[2];
rz(-0.52545737) q[3];
sx q[3];
rz(-0.94949978) q[3];
sx q[3];
rz(-0.8233101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7354108) q[2];
sx q[2];
rz(-1.7291131) q[2];
sx q[2];
rz(-2.0162876) q[2];
rz(-0.79648894) q[3];
sx q[3];
rz(-0.73553604) q[3];
sx q[3];
rz(-0.19851941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35561246) q[0];
sx q[0];
rz(-2.8481843) q[0];
sx q[0];
rz(-2.8134213) q[0];
rz(-2.7525821) q[1];
sx q[1];
rz(-1.5711454) q[1];
sx q[1];
rz(1.6860115) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7195014) q[0];
sx q[0];
rz(-1.4900643) q[0];
sx q[0];
rz(-2.9305365) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7057538) q[2];
sx q[2];
rz(-1.696387) q[2];
sx q[2];
rz(-1.3940982) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4715188) q[1];
sx q[1];
rz(-0.68738261) q[1];
sx q[1];
rz(1.9400466) q[1];
x q[2];
rz(1.9295272) q[3];
sx q[3];
rz(-1.609612) q[3];
sx q[3];
rz(-1.2286548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1652611) q[2];
sx q[2];
rz(-2.2165522) q[2];
sx q[2];
rz(0.46249214) q[2];
rz(-1.0235323) q[3];
sx q[3];
rz(-1.0547124) q[3];
sx q[3];
rz(0.83824497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8031215) q[0];
sx q[0];
rz(-1.7110889) q[0];
sx q[0];
rz(0.17369239) q[0];
rz(0.22131418) q[1];
sx q[1];
rz(-0.68562713) q[1];
sx q[1];
rz(1.3632704) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1817148) q[0];
sx q[0];
rz(-1.7522488) q[0];
sx q[0];
rz(-1.160301) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0603833) q[2];
sx q[2];
rz(-1.6575362) q[2];
sx q[2];
rz(-0.38151151) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2491571) q[1];
sx q[1];
rz(-1.1022864) q[1];
sx q[1];
rz(-2.4484642) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2962988) q[3];
sx q[3];
rz(-0.35652439) q[3];
sx q[3];
rz(-0.58024065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.26019105) q[2];
sx q[2];
rz(-1.9804201) q[2];
sx q[2];
rz(-1.6960404) q[2];
rz(0.77406231) q[3];
sx q[3];
rz(-2.4249707) q[3];
sx q[3];
rz(-1.4708446) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.562302) q[0];
sx q[0];
rz(-0.93920541) q[0];
sx q[0];
rz(0.24765177) q[0];
rz(2.472645) q[1];
sx q[1];
rz(-1.9525783) q[1];
sx q[1];
rz(-0.98562366) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0730608) q[0];
sx q[0];
rz(-1.4559696) q[0];
sx q[0];
rz(0.5590183) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.885576) q[2];
sx q[2];
rz(-1.0295964) q[2];
sx q[2];
rz(1.7371295) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0874153) q[1];
sx q[1];
rz(-1.2654516) q[1];
sx q[1];
rz(1.9806238) q[1];
rz(1.5721442) q[3];
sx q[3];
rz(-2.165757) q[3];
sx q[3];
rz(-1.6890989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.140427) q[2];
sx q[2];
rz(-0.75158921) q[2];
sx q[2];
rz(1.8482194) q[2];
rz(-1.9017396) q[3];
sx q[3];
rz(-2.445502) q[3];
sx q[3];
rz(2.4333439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86466113) q[0];
sx q[0];
rz(-1.4286574) q[0];
sx q[0];
rz(0.96555936) q[0];
rz(-2.7652265) q[1];
sx q[1];
rz(-1.8195567) q[1];
sx q[1];
rz(1.2300864) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9133643) q[0];
sx q[0];
rz(-1.3051864) q[0];
sx q[0];
rz(1.6227386) q[0];
rz(2.7509934) q[2];
sx q[2];
rz(-1.877344) q[2];
sx q[2];
rz(1.3424003) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9113637) q[1];
sx q[1];
rz(-1.5004284) q[1];
sx q[1];
rz(0.60649782) q[1];
rz(-pi) q[2];
rz(-3.0885124) q[3];
sx q[3];
rz(-2.3276668) q[3];
sx q[3];
rz(-0.9313213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.85912117) q[2];
sx q[2];
rz(-0.66994795) q[2];
sx q[2];
rz(-3.0100789) q[2];
rz(-0.48111835) q[3];
sx q[3];
rz(-0.93004623) q[3];
sx q[3];
rz(-2.6432945) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0216825) q[0];
sx q[0];
rz(-2.5639738) q[0];
sx q[0];
rz(2.6525894) q[0];
rz(-1.3267714) q[1];
sx q[1];
rz(-1.860447) q[1];
sx q[1];
rz(-2.9723523) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45825935) q[0];
sx q[0];
rz(-2.0067599) q[0];
sx q[0];
rz(-0.13157121) q[0];
rz(-pi) q[1];
rz(1.4812977) q[2];
sx q[2];
rz(-1.3907593) q[2];
sx q[2];
rz(-1.5370739) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0161017) q[1];
sx q[1];
rz(-0.91052848) q[1];
sx q[1];
rz(1.0339526) q[1];
rz(-pi) q[2];
rz(0.27948252) q[3];
sx q[3];
rz(-1.8516083) q[3];
sx q[3];
rz(-1.3407509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5883611) q[2];
sx q[2];
rz(-1.0855805) q[2];
sx q[2];
rz(-0.21151839) q[2];
rz(1.616098) q[3];
sx q[3];
rz(-1.0001405) q[3];
sx q[3];
rz(2.8741527) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78032988) q[0];
sx q[0];
rz(-1.7050671) q[0];
sx q[0];
rz(0.29062511) q[0];
rz(2.3296539) q[1];
sx q[1];
rz(-2.5135136) q[1];
sx q[1];
rz(1.5807349) q[1];
rz(1.590325) q[2];
sx q[2];
rz(-1.1371798) q[2];
sx q[2];
rz(-0.22424998) q[2];
rz(-0.91625253) q[3];
sx q[3];
rz(-0.8662681) q[3];
sx q[3];
rz(-2.7570799) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
