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
rz(-1.6260835) q[0];
sx q[0];
rz(-0.27684394) q[0];
sx q[0];
rz(-3.0247363) q[0];
rz(-0.46861831) q[1];
sx q[1];
rz(-2.7873971) q[1];
sx q[1];
rz(2.6503704) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6913585) q[0];
sx q[0];
rz(-1.629279) q[0];
sx q[0];
rz(-1.081715) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7560656) q[2];
sx q[2];
rz(-0.51947278) q[2];
sx q[2];
rz(-2.9105718) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3892606) q[1];
sx q[1];
rz(-1.6698784) q[1];
sx q[1];
rz(-2.5072751) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5003203) q[3];
sx q[3];
rz(-0.34409663) q[3];
sx q[3];
rz(-1.9650353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3789856) q[2];
sx q[2];
rz(-2.8540322) q[2];
sx q[2];
rz(-2.911705) q[2];
rz(-3.0843132) q[3];
sx q[3];
rz(-0.66904896) q[3];
sx q[3];
rz(2.2288286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7951935) q[0];
sx q[0];
rz(-0.38723463) q[0];
sx q[0];
rz(2.9246395) q[0];
rz(0.03288658) q[1];
sx q[1];
rz(-2.2078728) q[1];
sx q[1];
rz(0.65509534) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0467619) q[0];
sx q[0];
rz(-0.74709409) q[0];
sx q[0];
rz(1.8674891) q[0];
rz(-1.7704324) q[2];
sx q[2];
rz(-2.0642274) q[2];
sx q[2];
rz(0.1525998) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2456783) q[1];
sx q[1];
rz(-1.0686456) q[1];
sx q[1];
rz(-0.72523004) q[1];
rz(1.7294077) q[3];
sx q[3];
rz(-2.4008823) q[3];
sx q[3];
rz(-2.9155758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.81755) q[2];
sx q[2];
rz(-0.62411672) q[2];
sx q[2];
rz(-1.2980655) q[2];
rz(-1.921418) q[3];
sx q[3];
rz(-1.4379359) q[3];
sx q[3];
rz(-0.65742457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87797457) q[0];
sx q[0];
rz(-0.89732301) q[0];
sx q[0];
rz(-0.61400145) q[0];
rz(1.7738495) q[1];
sx q[1];
rz(-0.65679336) q[1];
sx q[1];
rz(2.6460389) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4456291) q[0];
sx q[0];
rz(-2.3167452) q[0];
sx q[0];
rz(2.4936409) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88286384) q[2];
sx q[2];
rz(-1.9830215) q[2];
sx q[2];
rz(-0.26160535) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.12901356) q[1];
sx q[1];
rz(-1.3316175) q[1];
sx q[1];
rz(-2.9496179) q[1];
rz(-pi) q[2];
rz(-0.63276799) q[3];
sx q[3];
rz(-0.55857165) q[3];
sx q[3];
rz(-2.4726923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.76501781) q[2];
sx q[2];
rz(-1.0708662) q[2];
sx q[2];
rz(0.59256727) q[2];
rz(2.9060034) q[3];
sx q[3];
rz(-2.3841136) q[3];
sx q[3];
rz(-0.45430115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6465103) q[0];
sx q[0];
rz(-2.3401234) q[0];
sx q[0];
rz(-3.1357646) q[0];
rz(0.5799154) q[1];
sx q[1];
rz(-0.35570759) q[1];
sx q[1];
rz(-0.52111202) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.482655) q[0];
sx q[0];
rz(-1.577943) q[0];
sx q[0];
rz(-0.004645017) q[0];
rz(-pi) q[1];
rz(0.96984152) q[2];
sx q[2];
rz(-1.5950987) q[2];
sx q[2];
rz(-1.5497335) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.5856984) q[1];
sx q[1];
rz(-2.0221077) q[1];
sx q[1];
rz(-1.1232308) q[1];
x q[2];
rz(0.31217869) q[3];
sx q[3];
rz(-1.494656) q[3];
sx q[3];
rz(-2.5654441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2516605) q[2];
sx q[2];
rz(-0.53091383) q[2];
sx q[2];
rz(2.8233675) q[2];
rz(0.71277726) q[3];
sx q[3];
rz(-2.8409499) q[3];
sx q[3];
rz(1.6596863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7707959) q[0];
sx q[0];
rz(-0.21368055) q[0];
sx q[0];
rz(-3.1053012) q[0];
rz(-2.6151784) q[1];
sx q[1];
rz(-1.4585835) q[1];
sx q[1];
rz(-0.2821736) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53647868) q[0];
sx q[0];
rz(-1.8465586) q[0];
sx q[0];
rz(-2.2173397) q[0];
rz(2.2065978) q[2];
sx q[2];
rz(-0.80342573) q[2];
sx q[2];
rz(-2.406081) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6101004) q[1];
sx q[1];
rz(-0.86750114) q[1];
sx q[1];
rz(-1.6628357) q[1];
rz(-pi) q[2];
rz(-0.75936486) q[3];
sx q[3];
rz(-0.96624422) q[3];
sx q[3];
rz(0.60455632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7669749) q[2];
sx q[2];
rz(-0.66156113) q[2];
sx q[2];
rz(-3.0401163) q[2];
rz(0.96595079) q[3];
sx q[3];
rz(-2.2209397) q[3];
sx q[3];
rz(3.0275893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58979708) q[0];
sx q[0];
rz(-2.9647201) q[0];
sx q[0];
rz(-1.9675323) q[0];
rz(-2.7575098) q[1];
sx q[1];
rz(-1.1840772) q[1];
sx q[1];
rz(-0.039965872) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46972457) q[0];
sx q[0];
rz(-1.7789808) q[0];
sx q[0];
rz(-1.8987659) q[0];
rz(-pi) q[1];
rz(-2.2117679) q[2];
sx q[2];
rz(-1.8397959) q[2];
sx q[2];
rz(-2.9421634) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74943042) q[1];
sx q[1];
rz(-2.431173) q[1];
sx q[1];
rz(2.602052) q[1];
x q[2];
rz(0.82620718) q[3];
sx q[3];
rz(-1.7306657) q[3];
sx q[3];
rz(-1.9217971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.93969718) q[2];
sx q[2];
rz(-2.5494718) q[2];
sx q[2];
rz(-2.3218018) q[2];
rz(0.78153265) q[3];
sx q[3];
rz(-2.235409) q[3];
sx q[3];
rz(-1.7614822) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49067295) q[0];
sx q[0];
rz(-1.0395721) q[0];
sx q[0];
rz(1.8248722) q[0];
rz(1.7735749) q[1];
sx q[1];
rz(-1.2317692) q[1];
sx q[1];
rz(-0.43693158) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9428944) q[0];
sx q[0];
rz(-1.6075396) q[0];
sx q[0];
rz(1.7612639) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75114366) q[2];
sx q[2];
rz(-0.1607543) q[2];
sx q[2];
rz(3.1243665) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.93561427) q[1];
sx q[1];
rz(-2.2956478) q[1];
sx q[1];
rz(1.0826675) q[1];
rz(-pi) q[2];
rz(-2.715492) q[3];
sx q[3];
rz(-1.0747194) q[3];
sx q[3];
rz(-0.29447281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.74066073) q[2];
sx q[2];
rz(-2.0621982) q[2];
sx q[2];
rz(-3.0069922) q[2];
rz(1.2234737) q[3];
sx q[3];
rz(-1.384602) q[3];
sx q[3];
rz(-1.1659762) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3400035) q[0];
sx q[0];
rz(-2.8497301) q[0];
sx q[0];
rz(-2.959751) q[0];
rz(0.34135154) q[1];
sx q[1];
rz(-2.0149442) q[1];
sx q[1];
rz(0.15886074) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6065935) q[0];
sx q[0];
rz(-2.6887021) q[0];
sx q[0];
rz(-0.15720982) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42271013) q[2];
sx q[2];
rz(-1.8313932) q[2];
sx q[2];
rz(0.24525951) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.73028683) q[1];
sx q[1];
rz(-1.3540233) q[1];
sx q[1];
rz(2.8764399) q[1];
rz(-2.8980632) q[3];
sx q[3];
rz(-0.56748828) q[3];
sx q[3];
rz(-1.6220038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2813256) q[2];
sx q[2];
rz(-2.211326) q[2];
sx q[2];
rz(-2.2567828) q[2];
rz(1.7224711) q[3];
sx q[3];
rz(-2.116674) q[3];
sx q[3];
rz(0.44986808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29174969) q[0];
sx q[0];
rz(-1.869864) q[0];
sx q[0];
rz(2.1480609) q[0];
rz(1.3604856) q[1];
sx q[1];
rz(-0.88881701) q[1];
sx q[1];
rz(2.5885168) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21453126) q[0];
sx q[0];
rz(-1.1916516) q[0];
sx q[0];
rz(-1.8168455) q[0];
x q[1];
rz(-0.07358609) q[2];
sx q[2];
rz(-2.4017757) q[2];
sx q[2];
rz(0.70895665) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83718761) q[1];
sx q[1];
rz(-0.68301979) q[1];
sx q[1];
rz(0.4153312) q[1];
rz(-1.2129799) q[3];
sx q[3];
rz(-2.4705437) q[3];
sx q[3];
rz(0.34085654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9370344) q[2];
sx q[2];
rz(-2.5356346) q[2];
sx q[2];
rz(-2.9045612) q[2];
rz(-1.5406746) q[3];
sx q[3];
rz(-2.377244) q[3];
sx q[3];
rz(-1.9717533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071526214) q[0];
sx q[0];
rz(-1.6306174) q[0];
sx q[0];
rz(-3.0226829) q[0];
rz(-0.43926829) q[1];
sx q[1];
rz(-1.2989137) q[1];
sx q[1];
rz(-0.063442245) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3812696) q[0];
sx q[0];
rz(-1.0372235) q[0];
sx q[0];
rz(-2.3593581) q[0];
x q[1];
rz(2.5173152) q[2];
sx q[2];
rz(-1.7221525) q[2];
sx q[2];
rz(0.69619149) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0234785) q[1];
sx q[1];
rz(-2.8148068) q[1];
sx q[1];
rz(-1.544373) q[1];
rz(2.2584553) q[3];
sx q[3];
rz(-2.4534907) q[3];
sx q[3];
rz(2.7795252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.1305609) q[2];
sx q[2];
rz(-2.06879) q[2];
sx q[2];
rz(1.6434742) q[2];
rz(0.9324075) q[3];
sx q[3];
rz(-1.1199718) q[3];
sx q[3];
rz(1.9278661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3907923) q[0];
sx q[0];
rz(-1.9600497) q[0];
sx q[0];
rz(2.0586769) q[0];
rz(-0.87286585) q[1];
sx q[1];
rz(-2.0582336) q[1];
sx q[1];
rz(2.2179926) q[1];
rz(0.93449705) q[2];
sx q[2];
rz(-2.1515982) q[2];
sx q[2];
rz(1.2610255) q[2];
rz(-0.57742248) q[3];
sx q[3];
rz(-1.631177) q[3];
sx q[3];
rz(-2.9384818) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
