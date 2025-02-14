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
rz(-1.5440829) q[0];
sx q[0];
rz(4.5151526) q[0];
sx q[0];
rz(8.406352) q[0];
rz(2.6234558) q[1];
sx q[1];
rz(5.5261373) q[1];
sx q[1];
rz(11.93461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.951183) q[0];
sx q[0];
rz(-2.8993239) q[0];
sx q[0];
rz(1.0231859) q[0];
rz(-1.6502981) q[2];
sx q[2];
rz(-0.30361807) q[2];
sx q[2];
rz(-2.7111766) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.5093928) q[1];
sx q[1];
rz(-0.94974564) q[1];
sx q[1];
rz(2.7953447) q[1];
x q[2];
rz(-1.8970892) q[3];
sx q[3];
rz(-1.4832116) q[3];
sx q[3];
rz(2.0503941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6112001) q[2];
sx q[2];
rz(-2.7860614) q[2];
sx q[2];
rz(1.4872888) q[2];
rz(-2.2330331) q[3];
sx q[3];
rz(-0.23659758) q[3];
sx q[3];
rz(-1.0049741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1210043) q[0];
sx q[0];
rz(-1.7089184) q[0];
sx q[0];
rz(-0.16227907) q[0];
rz(0.24457112) q[1];
sx q[1];
rz(-1.9385612) q[1];
sx q[1];
rz(-2.8499106) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26839089) q[0];
sx q[0];
rz(-1.8531728) q[0];
sx q[0];
rz(-1.3703521) q[0];
rz(-1.2873093) q[2];
sx q[2];
rz(-1.785673) q[2];
sx q[2];
rz(2.7708997) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6704935) q[1];
sx q[1];
rz(-2.3355977) q[1];
sx q[1];
rz(-0.32260311) q[1];
x q[2];
rz(0.79370768) q[3];
sx q[3];
rz(-2.5849403) q[3];
sx q[3];
rz(1.2964378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4633999) q[2];
sx q[2];
rz(-2.2726111) q[2];
sx q[2];
rz(-1.0758859) q[2];
rz(-1.2470657) q[3];
sx q[3];
rz(-1.5970634) q[3];
sx q[3];
rz(1.4472848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70188824) q[0];
sx q[0];
rz(-1.1264369) q[0];
sx q[0];
rz(-2.9578748) q[0];
rz(1.5048997) q[1];
sx q[1];
rz(-0.74438649) q[1];
sx q[1];
rz(0.16477975) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4713584) q[0];
sx q[0];
rz(-0.58588282) q[0];
sx q[0];
rz(-1.5878994) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8833227) q[2];
sx q[2];
rz(-2.222568) q[2];
sx q[2];
rz(1.5787909) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7565733) q[1];
sx q[1];
rz(-2.1655271) q[1];
sx q[1];
rz(-1.1043332) q[1];
rz(-pi) q[2];
rz(1.0298877) q[3];
sx q[3];
rz(-2.5383484) q[3];
sx q[3];
rz(-1.3689976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.1068716) q[2];
sx q[2];
rz(-1.160459) q[2];
sx q[2];
rz(-1.1064233) q[2];
rz(0.70118457) q[3];
sx q[3];
rz(-1.3283575) q[3];
sx q[3];
rz(0.4655233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.768854) q[0];
sx q[0];
rz(-2.0059858) q[0];
sx q[0];
rz(-0.94451529) q[0];
rz(-1.6230029) q[1];
sx q[1];
rz(-2.1606162) q[1];
sx q[1];
rz(1.6015582) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1683274) q[0];
sx q[0];
rz(-1.0078609) q[0];
sx q[0];
rz(-0.8416147) q[0];
x q[1];
rz(2.3643199) q[2];
sx q[2];
rz(-2.9748355) q[2];
sx q[2];
rz(0.47698944) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.67805144) q[1];
sx q[1];
rz(-1.1992595) q[1];
sx q[1];
rz(2.3080565) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.047961162) q[3];
sx q[3];
rz(-0.75100079) q[3];
sx q[3];
rz(-2.6329071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.010178415) q[2];
sx q[2];
rz(-0.62128908) q[2];
sx q[2];
rz(-0.16769257) q[2];
rz(0.0532648) q[3];
sx q[3];
rz(-2.0361418) q[3];
sx q[3];
rz(-1.3767327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5662956) q[0];
sx q[0];
rz(-3.0360041) q[0];
sx q[0];
rz(0.29933023) q[0];
rz(-2.7686367) q[1];
sx q[1];
rz(-1.8920218) q[1];
sx q[1];
rz(-1.4704871) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0051992) q[0];
sx q[0];
rz(-1.0233425) q[0];
sx q[0];
rz(-2.0776847) q[0];
x q[1];
rz(-1.4506571) q[2];
sx q[2];
rz(-1.1098301) q[2];
sx q[2];
rz(1.9751939) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.59488397) q[1];
sx q[1];
rz(-2.1404033) q[1];
sx q[1];
rz(0.96086603) q[1];
x q[2];
rz(-2.9657508) q[3];
sx q[3];
rz(-2.1851563) q[3];
sx q[3];
rz(-1.1509488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2028929) q[2];
sx q[2];
rz(-2.4788269) q[2];
sx q[2];
rz(-1.1104442) q[2];
rz(2.2796196) q[3];
sx q[3];
rz(-1.6317261) q[3];
sx q[3];
rz(-1.2601669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92457572) q[0];
sx q[0];
rz(-2.45455) q[0];
sx q[0];
rz(2.7365015) q[0];
rz(0.336126) q[1];
sx q[1];
rz(-1.4210217) q[1];
sx q[1];
rz(0.95132557) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2971056) q[0];
sx q[0];
rz(-1.1702303) q[0];
sx q[0];
rz(-1.2933008) q[0];
x q[1];
rz(-1.7826005) q[2];
sx q[2];
rz(-2.2957605) q[2];
sx q[2];
rz(-0.69378366) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3555043) q[1];
sx q[1];
rz(-2.5230683) q[1];
sx q[1];
rz(0.39822802) q[1];
rz(0.59744617) q[3];
sx q[3];
rz(-2.8077112) q[3];
sx q[3];
rz(-2.2528439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0508017) q[2];
sx q[2];
rz(-1.7369221) q[2];
sx q[2];
rz(-2.9108099) q[2];
rz(2.9595621) q[3];
sx q[3];
rz(-0.56806505) q[3];
sx q[3];
rz(0.52869421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.35293216) q[0];
sx q[0];
rz(-0.039529888) q[0];
sx q[0];
rz(-2.2094862) q[0];
rz(-2.508029) q[1];
sx q[1];
rz(-2.0220058) q[1];
sx q[1];
rz(-1.8036141) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.103823) q[0];
sx q[0];
rz(-1.6523721) q[0];
sx q[0];
rz(1.4647116) q[0];
rz(0.39057486) q[2];
sx q[2];
rz(-2.8468067) q[2];
sx q[2];
rz(0.24402555) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6465969) q[1];
sx q[1];
rz(-0.75237583) q[1];
sx q[1];
rz(2.1667891) q[1];
x q[2];
rz(-0.058480992) q[3];
sx q[3];
rz(-2.6425458) q[3];
sx q[3];
rz(-0.24802314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46447095) q[2];
sx q[2];
rz(-1.9147583) q[2];
sx q[2];
rz(1.9390437) q[2];
rz(-0.36457148) q[3];
sx q[3];
rz(-2.4489844) q[3];
sx q[3];
rz(-2.306849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4778336) q[0];
sx q[0];
rz(-0.57357016) q[0];
sx q[0];
rz(-0.17663503) q[0];
rz(0.44003507) q[1];
sx q[1];
rz(-1.1853848) q[1];
sx q[1];
rz(2.1766591) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7149219) q[0];
sx q[0];
rz(-1.3765125) q[0];
sx q[0];
rz(1.4358836) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72146031) q[2];
sx q[2];
rz(-2.4000492) q[2];
sx q[2];
rz(-2.3539345) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0256465) q[1];
sx q[1];
rz(-1.0315511) q[1];
sx q[1];
rz(0.42614062) q[1];
rz(-pi) q[2];
rz(0.075841622) q[3];
sx q[3];
rz(-1.8245398) q[3];
sx q[3];
rz(0.50204078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9289916) q[2];
sx q[2];
rz(-1.6182199) q[2];
sx q[2];
rz(-0.0099445899) q[2];
rz(-0.11659226) q[3];
sx q[3];
rz(-0.3392342) q[3];
sx q[3];
rz(-0.54178437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56124878) q[0];
sx q[0];
rz(-1.2224226) q[0];
sx q[0];
rz(-2.5196581) q[0];
rz(-1.4382582) q[1];
sx q[1];
rz(-2.56918) q[1];
sx q[1];
rz(-0.66351801) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7166876) q[0];
sx q[0];
rz(-0.7640673) q[0];
sx q[0];
rz(-0.31830799) q[0];
x q[1];
rz(-1.1691888) q[2];
sx q[2];
rz(-0.3065232) q[2];
sx q[2];
rz(-1.0108794) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.864526) q[1];
sx q[1];
rz(-0.73523318) q[1];
sx q[1];
rz(2.0551873) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28182192) q[3];
sx q[3];
rz(-1.661206) q[3];
sx q[3];
rz(-0.91109959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5180987) q[2];
sx q[2];
rz(-1.3891209) q[2];
sx q[2];
rz(-0.36250472) q[2];
rz(-0.47719657) q[3];
sx q[3];
rz(-1.0537182) q[3];
sx q[3];
rz(-1.1227192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057673205) q[0];
sx q[0];
rz(-2.4102983) q[0];
sx q[0];
rz(-2.2398563) q[0];
rz(2.4665191) q[1];
sx q[1];
rz(-2.156064) q[1];
sx q[1];
rz(0.14428446) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1146204) q[0];
sx q[0];
rz(-1.2204224) q[0];
sx q[0];
rz(0.087894126) q[0];
rz(-pi) q[1];
rz(2.921943) q[2];
sx q[2];
rz(-1.4712508) q[2];
sx q[2];
rz(1.2964028) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.54507885) q[1];
sx q[1];
rz(-2.1987913) q[1];
sx q[1];
rz(1.2596115) q[1];
rz(-pi) q[2];
rz(0.38548174) q[3];
sx q[3];
rz(-0.47244888) q[3];
sx q[3];
rz(2.2466898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5054063) q[2];
sx q[2];
rz(-1.9289086) q[2];
sx q[2];
rz(-0.20067659) q[2];
rz(1.1096654) q[3];
sx q[3];
rz(-2.6789013) q[3];
sx q[3];
rz(2.5698575) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5889482) q[0];
sx q[0];
rz(-0.968796) q[0];
sx q[0];
rz(2.0198685) q[0];
rz(-0.48925346) q[1];
sx q[1];
rz(-1.6901292) q[1];
sx q[1];
rz(2.1314175) q[1];
rz(-3.048754) q[2];
sx q[2];
rz(-1.1576817) q[2];
sx q[2];
rz(0.67410034) q[2];
rz(-0.44224593) q[3];
sx q[3];
rz(-1.5968235) q[3];
sx q[3];
rz(2.1376283) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
