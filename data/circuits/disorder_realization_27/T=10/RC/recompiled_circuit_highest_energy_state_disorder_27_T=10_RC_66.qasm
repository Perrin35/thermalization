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
rz(2.8994695) q[0];
sx q[0];
rz(-2.4162633) q[0];
sx q[0];
rz(2.2348833) q[0];
rz(-0.46696219) q[1];
sx q[1];
rz(4.0443647) q[1];
sx q[1];
rz(9.180896) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42729331) q[0];
sx q[0];
rz(-0.9511982) q[0];
sx q[0];
rz(-1.0788973) q[0];
rz(2.4588434) q[2];
sx q[2];
rz(-1.054316) q[2];
sx q[2];
rz(-0.7841332) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.697487) q[1];
sx q[1];
rz(-0.19679697) q[1];
sx q[1];
rz(-2.7641731) q[1];
x q[2];
rz(2.3238691) q[3];
sx q[3];
rz(-1.4151038) q[3];
sx q[3];
rz(-0.45785357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4563518) q[2];
sx q[2];
rz(-0.46619236) q[2];
sx q[2];
rz(-1.0249798) q[2];
rz(-3.0023365) q[3];
sx q[3];
rz(-1.9459629) q[3];
sx q[3];
rz(2.9774184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76992947) q[0];
sx q[0];
rz(-2.0491845) q[0];
sx q[0];
rz(1.9364233) q[0];
rz(-1.4681762) q[1];
sx q[1];
rz(-1.2638777) q[1];
sx q[1];
rz(-1.7211627) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.159585) q[0];
sx q[0];
rz(-0.8340584) q[0];
sx q[0];
rz(-2.9247051) q[0];
rz(-pi) q[1];
rz(-1.7718138) q[2];
sx q[2];
rz(-2.8716785) q[2];
sx q[2];
rz(2.6156099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0542595) q[1];
sx q[1];
rz(-1.9078322) q[1];
sx q[1];
rz(0.23566206) q[1];
rz(-pi) q[2];
rz(-1.0923493) q[3];
sx q[3];
rz(-0.88913267) q[3];
sx q[3];
rz(-2.9605856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.18985192) q[2];
sx q[2];
rz(-0.18988374) q[2];
sx q[2];
rz(-0.80113775) q[2];
rz(2.7028911) q[3];
sx q[3];
rz(-1.2480241) q[3];
sx q[3];
rz(-1.1130822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7640215) q[0];
sx q[0];
rz(-2.8442597) q[0];
sx q[0];
rz(1.1391033) q[0];
rz(2.8746919) q[1];
sx q[1];
rz(-1.2511988) q[1];
sx q[1];
rz(2.7762754) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034780033) q[0];
sx q[0];
rz(-1.0584373) q[0];
sx q[0];
rz(1.4573291) q[0];
x q[1];
rz(-2.2174805) q[2];
sx q[2];
rz(-2.036493) q[2];
sx q[2];
rz(0.52801029) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68054861) q[1];
sx q[1];
rz(-0.37255424) q[1];
sx q[1];
rz(1.4112605) q[1];
rz(-0.85478304) q[3];
sx q[3];
rz(-1.0589335) q[3];
sx q[3];
rz(1.8292242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.9109362) q[2];
sx q[2];
rz(-0.74030423) q[2];
sx q[2];
rz(-0.12348565) q[2];
rz(2.2423045) q[3];
sx q[3];
rz(-1.6920009) q[3];
sx q[3];
rz(1.2353157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7742291) q[0];
sx q[0];
rz(-1.4643359) q[0];
sx q[0];
rz(-0.88957077) q[0];
rz(2.3956237) q[1];
sx q[1];
rz(-1.6791226) q[1];
sx q[1];
rz(1.3847345) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9002514) q[0];
sx q[0];
rz(-1.6286193) q[0];
sx q[0];
rz(3.1067294) q[0];
rz(-pi) q[1];
rz(0.32082243) q[2];
sx q[2];
rz(-2.4745382) q[2];
sx q[2];
rz(-0.49122444) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3900323) q[1];
sx q[1];
rz(-0.96876013) q[1];
sx q[1];
rz(1.3436418) q[1];
rz(-pi) q[2];
rz(1.6024578) q[3];
sx q[3];
rz(-1.0370812) q[3];
sx q[3];
rz(-2.9759852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3978465) q[2];
sx q[2];
rz(-0.83796871) q[2];
sx q[2];
rz(-2.1261334) q[2];
rz(0.022653496) q[3];
sx q[3];
rz(-1.7706784) q[3];
sx q[3];
rz(3.0034351) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7669693) q[0];
sx q[0];
rz(-1.2767108) q[0];
sx q[0];
rz(-2.936777) q[0];
rz(-2.9264033) q[1];
sx q[1];
rz(-2.9209825) q[1];
sx q[1];
rz(-2.2664216) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29680934) q[0];
sx q[0];
rz(-1.267485) q[0];
sx q[0];
rz(-0.15547296) q[0];
rz(-pi) q[1];
rz(0.83568253) q[2];
sx q[2];
rz(-1.6854523) q[2];
sx q[2];
rz(-2.3029652) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2332747) q[1];
sx q[1];
rz(-1.868778) q[1];
sx q[1];
rz(1.3402935) q[1];
rz(0.93980333) q[3];
sx q[3];
rz(-1.6865589) q[3];
sx q[3];
rz(-0.94652688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.63336664) q[2];
sx q[2];
rz(-0.88819155) q[2];
sx q[2];
rz(1.4166895) q[2];
rz(2.21777) q[3];
sx q[3];
rz(-0.97771907) q[3];
sx q[3];
rz(-2.0098604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0210339) q[0];
sx q[0];
rz(-2.4036305) q[0];
sx q[0];
rz(2.293204) q[0];
rz(1.5658762) q[1];
sx q[1];
rz(-1.8137685) q[1];
sx q[1];
rz(0.94589344) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1740055) q[0];
sx q[0];
rz(-1.286309) q[0];
sx q[0];
rz(0.6993642) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81056548) q[2];
sx q[2];
rz(-1.422035) q[2];
sx q[2];
rz(-2.8523977) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0578445) q[1];
sx q[1];
rz(-1.9942772) q[1];
sx q[1];
rz(-2.4992056) q[1];
rz(-1.7686144) q[3];
sx q[3];
rz(-2.0536011) q[3];
sx q[3];
rz(1.5316666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9653603) q[2];
sx q[2];
rz(-1.2364028) q[2];
sx q[2];
rz(1.4005093) q[2];
rz(1.4817574) q[3];
sx q[3];
rz(-1.3503431) q[3];
sx q[3];
rz(-0.19937936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(1.4014715) q[0];
sx q[0];
rz(-2.6417702) q[0];
sx q[0];
rz(1.9299782) q[0];
rz(2.6648193) q[1];
sx q[1];
rz(-1.2766653) q[1];
sx q[1];
rz(1.5062987) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7619373) q[0];
sx q[0];
rz(-2.5003644) q[0];
sx q[0];
rz(-1.0851218) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25038211) q[2];
sx q[2];
rz(-1.2571475) q[2];
sx q[2];
rz(0.55082441) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9730649) q[1];
sx q[1];
rz(-2.2650411) q[1];
sx q[1];
rz(-1.7315052) q[1];
x q[2];
rz(-2.659306) q[3];
sx q[3];
rz(-1.3083754) q[3];
sx q[3];
rz(-0.69604128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3393042) q[2];
sx q[2];
rz(-1.5172639) q[2];
sx q[2];
rz(-0.18292546) q[2];
rz(-0.66169345) q[3];
sx q[3];
rz(-1.9293834) q[3];
sx q[3];
rz(2.4367512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45109192) q[0];
sx q[0];
rz(-1.4649614) q[0];
sx q[0];
rz(0.14725421) q[0];
rz(0.14689771) q[1];
sx q[1];
rz(-0.92120996) q[1];
sx q[1];
rz(-1.1796835) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058181949) q[0];
sx q[0];
rz(-1.9866052) q[0];
sx q[0];
rz(-1.9394298) q[0];
rz(-pi) q[1];
rz(-0.4381035) q[2];
sx q[2];
rz(-2.0102276) q[2];
sx q[2];
rz(-0.319744) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.50936401) q[1];
sx q[1];
rz(-1.6809789) q[1];
sx q[1];
rz(-2.5446507) q[1];
x q[2];
rz(-0.75408077) q[3];
sx q[3];
rz(-1.4907661) q[3];
sx q[3];
rz(-2.4802993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94072056) q[2];
sx q[2];
rz(-2.2042037) q[2];
sx q[2];
rz(-2.3737523) q[2];
rz(-2.614894) q[3];
sx q[3];
rz(-2.0898762) q[3];
sx q[3];
rz(0.55829486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4621157) q[0];
sx q[0];
rz(-2.9556892) q[0];
sx q[0];
rz(2.0515077) q[0];
rz(-1.3500301) q[1];
sx q[1];
rz(-1.7410024) q[1];
sx q[1];
rz(-2.7166691) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7658378) q[0];
sx q[0];
rz(-1.100044) q[0];
sx q[0];
rz(0.063675675) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64984729) q[2];
sx q[2];
rz(-1.6260626) q[2];
sx q[2];
rz(-1.2853607) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.35985816) q[1];
sx q[1];
rz(-1.1175851) q[1];
sx q[1];
rz(1.9537326) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70885371) q[3];
sx q[3];
rz(-1.1579305) q[3];
sx q[3];
rz(2.9065437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2316042) q[2];
sx q[2];
rz(-0.51077545) q[2];
sx q[2];
rz(-2.6166022) q[2];
rz(1.7867583) q[3];
sx q[3];
rz(-0.89559186) q[3];
sx q[3];
rz(1.7790599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7056284) q[0];
sx q[0];
rz(-2.4810677) q[0];
sx q[0];
rz(0.10667644) q[0];
rz(1.0542997) q[1];
sx q[1];
rz(-1.432447) q[1];
sx q[1];
rz(-1.7342825) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5085691) q[0];
sx q[0];
rz(-2.3481524) q[0];
sx q[0];
rz(2.2587848) q[0];
rz(-pi) q[1];
rz(2.5564747) q[2];
sx q[2];
rz(-2.3986534) q[2];
sx q[2];
rz(-1.20005) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5104502) q[1];
sx q[1];
rz(-1.4179953) q[1];
sx q[1];
rz(1.3451206) q[1];
rz(-0.75026433) q[3];
sx q[3];
rz(-1.0862873) q[3];
sx q[3];
rz(1.110525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6537031) q[2];
sx q[2];
rz(-2.6370878) q[2];
sx q[2];
rz(0.27296909) q[2];
rz(-2.2702787) q[3];
sx q[3];
rz(-1.4600236) q[3];
sx q[3];
rz(3.0070846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3684261) q[0];
sx q[0];
rz(-1.9482524) q[0];
sx q[0];
rz(0.85309749) q[0];
rz(0.38289616) q[1];
sx q[1];
rz(-1.2877512) q[1];
sx q[1];
rz(3.126694) q[1];
rz(0.64528633) q[2];
sx q[2];
rz(-2.4748079) q[2];
sx q[2];
rz(2.9590068) q[2];
rz(-0.79991585) q[3];
sx q[3];
rz(-1.6408995) q[3];
sx q[3];
rz(-1.2560918) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
