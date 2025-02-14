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
rz(-2.146848) q[0];
sx q[0];
rz(-0.58726197) q[0];
sx q[0];
rz(2.5160312) q[0];
rz(0.62043959) q[1];
sx q[1];
rz(4.1320463) q[1];
sx q[1];
rz(10.208701) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6990358) q[0];
sx q[0];
rz(-0.43414206) q[0];
sx q[0];
rz(3.1331557) q[0];
rz(-pi) q[1];
rz(-0.13909362) q[2];
sx q[2];
rz(-1.6533829) q[2];
sx q[2];
rz(2.6001584) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82134889) q[1];
sx q[1];
rz(-0.87961266) q[1];
sx q[1];
rz(2.3752179) q[1];
rz(-pi) q[2];
rz(-1.6585856) q[3];
sx q[3];
rz(-0.79500073) q[3];
sx q[3];
rz(2.3203497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94773942) q[2];
sx q[2];
rz(-2.096602) q[2];
sx q[2];
rz(-2.4742773) q[2];
rz(0.31467485) q[3];
sx q[3];
rz(-2.5421725) q[3];
sx q[3];
rz(-0.92710322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.2941403) q[0];
sx q[0];
rz(-0.84646928) q[0];
sx q[0];
rz(2.8417929) q[0];
rz(-1.0046129) q[1];
sx q[1];
rz(-0.71669465) q[1];
sx q[1];
rz(-2.2915548) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3687821) q[0];
sx q[0];
rz(-1.6666935) q[0];
sx q[0];
rz(2.1762288) q[0];
rz(-pi) q[1];
rz(-2.9710567) q[2];
sx q[2];
rz(-1.6493622) q[2];
sx q[2];
rz(0.54135347) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5932754) q[1];
sx q[1];
rz(-1.2139582) q[1];
sx q[1];
rz(0.22290454) q[1];
x q[2];
rz(-1.5845234) q[3];
sx q[3];
rz(-1.9753595) q[3];
sx q[3];
rz(1.8006397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5689508) q[2];
sx q[2];
rz(-2.6093542) q[2];
sx q[2];
rz(0.66298318) q[2];
rz(0.7524544) q[3];
sx q[3];
rz(-1.821725) q[3];
sx q[3];
rz(-2.9451356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74333423) q[0];
sx q[0];
rz(-2.7404116) q[0];
sx q[0];
rz(2.4617526) q[0];
rz(-2.4196449) q[1];
sx q[1];
rz(-0.55689055) q[1];
sx q[1];
rz(0.78071761) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3975383) q[0];
sx q[0];
rz(-0.40073943) q[0];
sx q[0];
rz(-1.5512054) q[0];
rz(-pi) q[1];
rz(-0.0038006101) q[2];
sx q[2];
rz(-1.4846621) q[2];
sx q[2];
rz(-2.0250367) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.225395) q[1];
sx q[1];
rz(-0.94630264) q[1];
sx q[1];
rz(-1.1033113) q[1];
rz(-pi) q[2];
rz(-2.3060477) q[3];
sx q[3];
rz(-1.4073331) q[3];
sx q[3];
rz(-0.36629656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9146933) q[2];
sx q[2];
rz(-1.9136027) q[2];
sx q[2];
rz(2.4719888) q[2];
rz(-2.2412444) q[3];
sx q[3];
rz(-1.0122296) q[3];
sx q[3];
rz(-2.3646234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25903073) q[0];
sx q[0];
rz(-0.43864033) q[0];
sx q[0];
rz(-0.90748179) q[0];
rz(-0.59440815) q[1];
sx q[1];
rz(-0.92955697) q[1];
sx q[1];
rz(1.1297191) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66624712) q[0];
sx q[0];
rz(-2.1985558) q[0];
sx q[0];
rz(0.23601377) q[0];
x q[1];
rz(1.8215034) q[2];
sx q[2];
rz(-0.54566979) q[2];
sx q[2];
rz(0.16599338) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.10953021) q[1];
sx q[1];
rz(-2.1550165) q[1];
sx q[1];
rz(-1.8971127) q[1];
rz(-2.8606671) q[3];
sx q[3];
rz(-1.2884029) q[3];
sx q[3];
rz(-0.97746935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0323459) q[2];
sx q[2];
rz(-1.3999908) q[2];
sx q[2];
rz(-2.4692811) q[2];
rz(0.0191056) q[3];
sx q[3];
rz(-0.34135434) q[3];
sx q[3];
rz(-3.0066971) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8371209) q[0];
sx q[0];
rz(-0.79455513) q[0];
sx q[0];
rz(-0.16803148) q[0];
rz(1.2270323) q[1];
sx q[1];
rz(-1.2201759) q[1];
sx q[1];
rz(2.8438445) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15140238) q[0];
sx q[0];
rz(-2.08648) q[0];
sx q[0];
rz(-0.97008743) q[0];
rz(-pi) q[1];
rz(2.3565294) q[2];
sx q[2];
rz(-2.3975971) q[2];
sx q[2];
rz(0.97563484) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9492901) q[1];
sx q[1];
rz(-1.5459035) q[1];
sx q[1];
rz(0.20583238) q[1];
x q[2];
rz(-1.566542) q[3];
sx q[3];
rz(-2.9828072) q[3];
sx q[3];
rz(0.53200737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1905404) q[2];
sx q[2];
rz(-2.0420044) q[2];
sx q[2];
rz(2.5229048) q[2];
rz(0.80858532) q[3];
sx q[3];
rz(-0.4777258) q[3];
sx q[3];
rz(1.3396858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65281463) q[0];
sx q[0];
rz(-1.5050911) q[0];
sx q[0];
rz(-2.6763647) q[0];
rz(0.48509994) q[1];
sx q[1];
rz(-0.41050375) q[1];
sx q[1];
rz(3.0533275) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3028835) q[0];
sx q[0];
rz(-2.0219936) q[0];
sx q[0];
rz(-2.124435) q[0];
rz(-3.1287304) q[2];
sx q[2];
rz(-2.2996827) q[2];
sx q[2];
rz(0.60573214) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3151206) q[1];
sx q[1];
rz(-1.406028) q[1];
sx q[1];
rz(2.3830448) q[1];
x q[2];
rz(-2.8537684) q[3];
sx q[3];
rz(-1.6555641) q[3];
sx q[3];
rz(2.1957977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7469067) q[2];
sx q[2];
rz(-2.2862819) q[2];
sx q[2];
rz(2.6439903) q[2];
rz(-0.48007128) q[3];
sx q[3];
rz(-2.2205455) q[3];
sx q[3];
rz(0.068643071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8467872) q[0];
sx q[0];
rz(-2.6755896) q[0];
sx q[0];
rz(-1.4303327) q[0];
rz(-1.3149186) q[1];
sx q[1];
rz(-2.7480835) q[1];
sx q[1];
rz(1.5865145) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1436227) q[0];
sx q[0];
rz(-0.41150948) q[0];
sx q[0];
rz(-1.8726148) q[0];
x q[1];
rz(2.5539923) q[2];
sx q[2];
rz(-0.94600224) q[2];
sx q[2];
rz(1.3252362) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79826971) q[1];
sx q[1];
rz(-0.99162356) q[1];
sx q[1];
rz(3.0208241) q[1];
rz(-pi) q[2];
rz(-1.6756546) q[3];
sx q[3];
rz(-0.63408454) q[3];
sx q[3];
rz(0.88632562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2985208) q[2];
sx q[2];
rz(-2.4303747) q[2];
sx q[2];
rz(-0.92740518) q[2];
rz(1.983042) q[3];
sx q[3];
rz(-0.68803334) q[3];
sx q[3];
rz(-3.0454175) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4226828) q[0];
sx q[0];
rz(-0.89253187) q[0];
sx q[0];
rz(0.476015) q[0];
rz(-0.99126518) q[1];
sx q[1];
rz(-0.74949336) q[1];
sx q[1];
rz(3.057726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059444025) q[0];
sx q[0];
rz(-2.5184439) q[0];
sx q[0];
rz(2.2150458) q[0];
x q[1];
rz(0.65292439) q[2];
sx q[2];
rz(-2.5819375) q[2];
sx q[2];
rz(-1.6728354) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.36271154) q[1];
sx q[1];
rz(-1.485752) q[1];
sx q[1];
rz(-1.9910732) q[1];
x q[2];
rz(-1.2060542) q[3];
sx q[3];
rz(-1.6514773) q[3];
sx q[3];
rz(-1.9239359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2241761) q[2];
sx q[2];
rz(-2.8755964) q[2];
sx q[2];
rz(2.4931397) q[2];
rz(2.2367541) q[3];
sx q[3];
rz(-1.6831393) q[3];
sx q[3];
rz(0.31074935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7998578) q[0];
sx q[0];
rz(-0.71001995) q[0];
sx q[0];
rz(-2.572686) q[0];
rz(-2.3078602) q[1];
sx q[1];
rz(-1.2833475) q[1];
sx q[1];
rz(2.1368829) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2434655) q[0];
sx q[0];
rz(-1.3782189) q[0];
sx q[0];
rz(2.5054727) q[0];
rz(1.1962074) q[2];
sx q[2];
rz(-0.56972144) q[2];
sx q[2];
rz(1.8853055) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8160287) q[1];
sx q[1];
rz(-0.70239151) q[1];
sx q[1];
rz(-1.4588474) q[1];
x q[2];
rz(-1.6528204) q[3];
sx q[3];
rz(-0.98482212) q[3];
sx q[3];
rz(0.88912933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7052762) q[2];
sx q[2];
rz(-2.269561) q[2];
sx q[2];
rz(2.7936068) q[2];
rz(-0.82585382) q[3];
sx q[3];
rz(-0.21819849) q[3];
sx q[3];
rz(-2.5364449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0843049) q[0];
sx q[0];
rz(-2.075752) q[0];
sx q[0];
rz(-0.2070981) q[0];
rz(0.53889489) q[1];
sx q[1];
rz(-2.5205044) q[1];
sx q[1];
rz(-2.5394687) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3823366) q[0];
sx q[0];
rz(-0.70092541) q[0];
sx q[0];
rz(0.78030326) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5670577) q[2];
sx q[2];
rz(-2.4705187) q[2];
sx q[2];
rz(0.35071638) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6844897) q[1];
sx q[1];
rz(-0.41375638) q[1];
sx q[1];
rz(1.3828418) q[1];
rz(-pi) q[2];
rz(-2.7902428) q[3];
sx q[3];
rz(-0.77882871) q[3];
sx q[3];
rz(0.035127775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3567317) q[2];
sx q[2];
rz(-0.95180231) q[2];
sx q[2];
rz(1.3447364) q[2];
rz(-0.52062672) q[3];
sx q[3];
rz(-0.27025637) q[3];
sx q[3];
rz(2.3394913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6945334) q[0];
sx q[0];
rz(-0.72493989) q[0];
sx q[0];
rz(-1.3865393) q[0];
rz(-2.3583892) q[1];
sx q[1];
rz(-1.422311) q[1];
sx q[1];
rz(1.5625988) q[1];
rz(1.8405909) q[2];
sx q[2];
rz(-1.9023583) q[2];
sx q[2];
rz(-2.6826774) q[2];
rz(1.5824885) q[3];
sx q[3];
rz(-1.3306918) q[3];
sx q[3];
rz(-0.065879475) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
