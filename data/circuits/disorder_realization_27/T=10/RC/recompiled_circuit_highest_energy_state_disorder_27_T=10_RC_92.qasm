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
rz(-0.24212317) q[0];
sx q[0];
rz(5.557856) q[0];
sx q[0];
rz(10.331487) q[0];
rz(-0.46696219) q[1];
sx q[1];
rz(-2.2388206) q[1];
sx q[1];
rz(-0.24388193) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31793247) q[0];
sx q[0];
rz(-0.77046227) q[0];
sx q[0];
rz(0.58485683) q[0];
x q[1];
rz(2.4087231) q[2];
sx q[2];
rz(-2.3112765) q[2];
sx q[2];
rz(2.900659) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.49748028) q[1];
sx q[1];
rz(-1.4986769) q[1];
sx q[1];
rz(0.18326541) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3238691) q[3];
sx q[3];
rz(-1.7264888) q[3];
sx q[3];
rz(-0.45785357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.68524086) q[2];
sx q[2];
rz(-0.46619236) q[2];
sx q[2];
rz(-2.1166128) q[2];
rz(0.13925615) q[3];
sx q[3];
rz(-1.9459629) q[3];
sx q[3];
rz(2.9774184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.3716632) q[0];
sx q[0];
rz(-1.0924082) q[0];
sx q[0];
rz(-1.2051693) q[0];
rz(-1.6734164) q[1];
sx q[1];
rz(-1.877715) q[1];
sx q[1];
rz(-1.7211627) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.842671) q[0];
sx q[0];
rz(-2.3793829) q[0];
sx q[0];
rz(1.803714) q[0];
rz(1.3697789) q[2];
sx q[2];
rz(-0.2699142) q[2];
sx q[2];
rz(-2.6156099) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.08733315) q[1];
sx q[1];
rz(-1.9078322) q[1];
sx q[1];
rz(2.9059306) q[1];
x q[2];
rz(-2.0492433) q[3];
sx q[3];
rz(-0.88913267) q[3];
sx q[3];
rz(2.9605856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9517407) q[2];
sx q[2];
rz(-2.9517089) q[2];
sx q[2];
rz(2.3404549) q[2];
rz(0.43870157) q[3];
sx q[3];
rz(-1.2480241) q[3];
sx q[3];
rz(-2.0285105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7640215) q[0];
sx q[0];
rz(-2.8442597) q[0];
sx q[0];
rz(1.1391033) q[0];
rz(-0.26690075) q[1];
sx q[1];
rz(-1.8903939) q[1];
sx q[1];
rz(0.36531726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6613839) q[0];
sx q[0];
rz(-1.669642) q[0];
sx q[0];
rz(0.5151202) q[0];
rz(-0.56198012) q[2];
sx q[2];
rz(-1.0022853) q[2];
sx q[2];
rz(-1.3696483) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.74148948) q[1];
sx q[1];
rz(-1.5129397) q[1];
sx q[1];
rz(-1.9390381) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5015058) q[3];
sx q[3];
rz(-0.96154204) q[3];
sx q[3];
rz(-2.9971788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2306564) q[2];
sx q[2];
rz(-0.74030423) q[2];
sx q[2];
rz(3.018107) q[2];
rz(-2.2423045) q[3];
sx q[3];
rz(-1.6920009) q[3];
sx q[3];
rz(-1.2353157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3673636) q[0];
sx q[0];
rz(-1.6772567) q[0];
sx q[0];
rz(2.2520219) q[0];
rz(-0.74596897) q[1];
sx q[1];
rz(-1.6791226) q[1];
sx q[1];
rz(1.3847345) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.810122) q[0];
sx q[0];
rz(-1.6056013) q[0];
sx q[0];
rz(-1.5129382) q[0];
rz(-pi) q[1];
rz(1.3273932) q[2];
sx q[2];
rz(-0.94331532) q[2];
sx q[2];
rz(0.091077591) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1390723) q[1];
sx q[1];
rz(-2.5031282) q[1];
sx q[1];
rz(-2.8248766) q[1];
rz(-pi) q[2];
rz(-1.5391348) q[3];
sx q[3];
rz(-1.0370812) q[3];
sx q[3];
rz(-2.9759852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3978465) q[2];
sx q[2];
rz(-0.83796871) q[2];
sx q[2];
rz(-1.0154593) q[2];
rz(0.022653496) q[3];
sx q[3];
rz(-1.3709143) q[3];
sx q[3];
rz(-3.0034351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7669693) q[0];
sx q[0];
rz(-1.8648819) q[0];
sx q[0];
rz(2.936777) q[0];
rz(-2.9264033) q[1];
sx q[1];
rz(-2.9209825) q[1];
sx q[1];
rz(0.87517103) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29680934) q[0];
sx q[0];
rz(-1.8741076) q[0];
sx q[0];
rz(-0.15547296) q[0];
x q[1];
rz(2.9875675) q[2];
sx q[2];
rz(-0.8416033) q[2];
sx q[2];
rz(2.306349) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.73126792) q[1];
sx q[1];
rz(-1.7909697) q[1];
sx q[1];
rz(-0.30558821) q[1];
x q[2];
rz(-2.9985626) q[3];
sx q[3];
rz(-2.1969079) q[3];
sx q[3];
rz(0.54007733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.508226) q[2];
sx q[2];
rz(-2.2534011) q[2];
sx q[2];
rz(-1.7249031) q[2];
rz(-0.9238227) q[3];
sx q[3];
rz(-2.1638736) q[3];
sx q[3];
rz(2.0098604) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0210339) q[0];
sx q[0];
rz(-0.73796213) q[0];
sx q[0];
rz(-2.293204) q[0];
rz(-1.5658762) q[1];
sx q[1];
rz(-1.8137685) q[1];
sx q[1];
rz(-0.94589344) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28117546) q[0];
sx q[0];
rz(-0.74587727) q[0];
sx q[0];
rz(-2.7152048) q[0];
rz(-pi) q[1];
rz(-1.7849493) q[2];
sx q[2];
rz(-0.77177599) q[2];
sx q[2];
rz(1.1269778) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1887363) q[1];
sx q[1];
rz(-2.1486001) q[1];
sx q[1];
rz(2.0835447) q[1];
rz(-pi) q[2];
rz(0.49094851) q[3];
sx q[3];
rz(-1.3958389) q[3];
sx q[3];
rz(-3.0096731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1762323) q[2];
sx q[2];
rz(-1.2364028) q[2];
sx q[2];
rz(-1.7410834) q[2];
rz(-1.4817574) q[3];
sx q[3];
rz(-1.7912495) q[3];
sx q[3];
rz(2.9422133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7401212) q[0];
sx q[0];
rz(-0.49982247) q[0];
sx q[0];
rz(-1.9299782) q[0];
rz(2.6648193) q[1];
sx q[1];
rz(-1.8649273) q[1];
sx q[1];
rz(1.635294) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20902987) q[0];
sx q[0];
rz(-1.8537921) q[0];
sx q[0];
rz(2.1542633) q[0];
rz(1.2477307) q[2];
sx q[2];
rz(-1.808721) q[2];
sx q[2];
rz(1.0987154) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9730649) q[1];
sx q[1];
rz(-0.87655156) q[1];
sx q[1];
rz(-1.7315052) q[1];
rz(-1.8651857) q[3];
sx q[3];
rz(-2.0352425) q[3];
sx q[3];
rz(2.4018231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8022884) q[2];
sx q[2];
rz(-1.5172639) q[2];
sx q[2];
rz(0.18292546) q[2];
rz(-2.4798992) q[3];
sx q[3];
rz(-1.9293834) q[3];
sx q[3];
rz(-2.4367512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45109192) q[0];
sx q[0];
rz(-1.6766312) q[0];
sx q[0];
rz(-2.9943384) q[0];
rz(2.9946949) q[1];
sx q[1];
rz(-0.92120996) q[1];
sx q[1];
rz(1.1796835) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0834107) q[0];
sx q[0];
rz(-1.9866052) q[0];
sx q[0];
rz(-1.2021628) q[0];
rz(-pi) q[1];
rz(-1.0919763) q[2];
sx q[2];
rz(-1.1767595) q[2];
sx q[2];
rz(2.0872598) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6322286) q[1];
sx q[1];
rz(-1.6809789) q[1];
sx q[1];
rz(0.59694196) q[1];
x q[2];
rz(0.75408077) q[3];
sx q[3];
rz(-1.4907661) q[3];
sx q[3];
rz(-0.66129337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.94072056) q[2];
sx q[2];
rz(-0.93738896) q[2];
sx q[2];
rz(2.3737523) q[2];
rz(-0.52669865) q[3];
sx q[3];
rz(-1.0517164) q[3];
sx q[3];
rz(-2.5832978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4621157) q[0];
sx q[0];
rz(-0.18590346) q[0];
sx q[0];
rz(-2.0515077) q[0];
rz(1.7915626) q[1];
sx q[1];
rz(-1.4005902) q[1];
sx q[1];
rz(2.7166691) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23608828) q[0];
sx q[0];
rz(-0.47471913) q[0];
sx q[0];
rz(1.4464054) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4917454) q[2];
sx q[2];
rz(-1.6260626) q[2];
sx q[2];
rz(1.2853607) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7577014) q[1];
sx q[1];
rz(-0.5846068) q[1];
sx q[1];
rz(-0.65442185) q[1];
x q[2];
rz(0.70885371) q[3];
sx q[3];
rz(-1.1579305) q[3];
sx q[3];
rz(2.9065437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9099884) q[2];
sx q[2];
rz(-0.51077545) q[2];
sx q[2];
rz(-0.52499047) q[2];
rz(-1.3548343) q[3];
sx q[3];
rz(-2.2460008) q[3];
sx q[3];
rz(-1.7790599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43596426) q[0];
sx q[0];
rz(-0.66052496) q[0];
sx q[0];
rz(-3.0349162) q[0];
rz(2.0872929) q[1];
sx q[1];
rz(-1.7091457) q[1];
sx q[1];
rz(-1.7342825) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76864734) q[0];
sx q[0];
rz(-2.1539186) q[0];
sx q[0];
rz(2.5685422) q[0];
rz(-pi) q[1];
rz(0.65347341) q[2];
sx q[2];
rz(-1.9536886) q[2];
sx q[2];
rz(0.083195638) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.095276621) q[1];
sx q[1];
rz(-1.3477948) q[1];
sx q[1];
rz(2.9848802) q[1];
rz(-pi) q[2];
rz(-2.484177) q[3];
sx q[3];
rz(-0.86674009) q[3];
sx q[3];
rz(0.92368607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4878896) q[2];
sx q[2];
rz(-0.50450486) q[2];
sx q[2];
rz(-0.27296909) q[2];
rz(2.2702787) q[3];
sx q[3];
rz(-1.4600236) q[3];
sx q[3];
rz(0.13450809) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7731666) q[0];
sx q[0];
rz(-1.9482524) q[0];
sx q[0];
rz(0.85309749) q[0];
rz(0.38289616) q[1];
sx q[1];
rz(-1.2877512) q[1];
sx q[1];
rz(3.126694) q[1];
rz(-0.64528633) q[2];
sx q[2];
rz(-0.66678478) q[2];
sx q[2];
rz(-0.18258584) q[2];
rz(-2.3416768) q[3];
sx q[3];
rz(-1.5006931) q[3];
sx q[3];
rz(1.8855008) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
