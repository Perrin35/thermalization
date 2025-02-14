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
rz(-0.023802726) q[0];
sx q[0];
rz(4.2476141) q[0];
sx q[0];
rz(10.201693) q[0];
rz(1.39224) q[1];
sx q[1];
rz(-1.3148146) q[1];
sx q[1];
rz(2.1652752) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2433246) q[0];
sx q[0];
rz(-1.7043566) q[0];
sx q[0];
rz(-1.1569174) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56371477) q[2];
sx q[2];
rz(-1.5140805) q[2];
sx q[2];
rz(-0.67981718) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3562505) q[1];
sx q[1];
rz(-1.9948729) q[1];
sx q[1];
rz(2.8431358) q[1];
rz(-pi) q[2];
rz(-0.99774811) q[3];
sx q[3];
rz(-1.5978509) q[3];
sx q[3];
rz(1.388035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6355847) q[2];
sx q[2];
rz(-3.0878461) q[2];
sx q[2];
rz(-0.40679833) q[2];
rz(-2.9721416) q[3];
sx q[3];
rz(-2.6120766) q[3];
sx q[3];
rz(2.0690401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2605543) q[0];
sx q[0];
rz(-0.23242234) q[0];
sx q[0];
rz(-0.0037923092) q[0];
rz(-3.0637528) q[1];
sx q[1];
rz(-0.65996116) q[1];
sx q[1];
rz(-0.30581623) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90214848) q[0];
sx q[0];
rz(-1.7551219) q[0];
sx q[0];
rz(-1.0499766) q[0];
rz(1.7051093) q[2];
sx q[2];
rz(-2.3181653) q[2];
sx q[2];
rz(3.1253377) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21286035) q[1];
sx q[1];
rz(-1.9668192) q[1];
sx q[1];
rz(-2.3235882) q[1];
rz(-3.0806957) q[3];
sx q[3];
rz(-1.506926) q[3];
sx q[3];
rz(-2.8012432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1514312) q[2];
sx q[2];
rz(-1.3913245) q[2];
sx q[2];
rz(-2.4988417) q[2];
rz(-0.07240545) q[3];
sx q[3];
rz(-1.0591155) q[3];
sx q[3];
rz(1.7806627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0533503) q[0];
sx q[0];
rz(-2.998816) q[0];
sx q[0];
rz(-2.9852168) q[0];
rz(-0.0414255) q[1];
sx q[1];
rz(-2.5138469) q[1];
sx q[1];
rz(1.5511537) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0254733) q[0];
sx q[0];
rz(-2.1120694) q[0];
sx q[0];
rz(-1.3295637) q[0];
rz(-pi) q[1];
rz(-0.65608187) q[2];
sx q[2];
rz(-1.6245884) q[2];
sx q[2];
rz(-2.057586) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9030021) q[1];
sx q[1];
rz(-0.5410453) q[1];
sx q[1];
rz(0.21958406) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4778792) q[3];
sx q[3];
rz(-2.1557689) q[3];
sx q[3];
rz(0.044540964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7330043) q[2];
sx q[2];
rz(-2.6027347) q[2];
sx q[2];
rz(0.020922529) q[2];
rz(2.9544592) q[3];
sx q[3];
rz(-0.20550607) q[3];
sx q[3];
rz(-3.0278897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1784096) q[0];
sx q[0];
rz(-2.8091176) q[0];
sx q[0];
rz(-2.7030429) q[0];
rz(1.6167538) q[1];
sx q[1];
rz(-2.8068145) q[1];
sx q[1];
rz(2.8964892) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1634451) q[0];
sx q[0];
rz(-2.4171099) q[0];
sx q[0];
rz(-1.4446769) q[0];
rz(-pi) q[1];
rz(2.7736362) q[2];
sx q[2];
rz(-1.7178255) q[2];
sx q[2];
rz(1.8859175) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72232258) q[1];
sx q[1];
rz(-0.55907226) q[1];
sx q[1];
rz(1.30617) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.887759) q[3];
sx q[3];
rz(-0.84268236) q[3];
sx q[3];
rz(-1.1945981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.601292) q[2];
sx q[2];
rz(-2.7056594) q[2];
sx q[2];
rz(0.33176804) q[2];
rz(-0.48745421) q[3];
sx q[3];
rz(-2.09477) q[3];
sx q[3];
rz(-2.148518) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29193923) q[0];
sx q[0];
rz(-1.6878457) q[0];
sx q[0];
rz(2.3680903) q[0];
rz(1.1812814) q[1];
sx q[1];
rz(-0.14207323) q[1];
sx q[1];
rz(-1.389651) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4163602) q[0];
sx q[0];
rz(-1.9342039) q[0];
sx q[0];
rz(-2.711722) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85948617) q[2];
sx q[2];
rz(-2.0469249) q[2];
sx q[2];
rz(1.7993594) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70848318) q[1];
sx q[1];
rz(-2.6603087) q[1];
sx q[1];
rz(2.9631056) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64719154) q[3];
sx q[3];
rz(-2.8464937) q[3];
sx q[3];
rz(1.5403252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2191849) q[2];
sx q[2];
rz(-1.9318523) q[2];
sx q[2];
rz(-0.47214559) q[2];
rz(1.2989429) q[3];
sx q[3];
rz(-1.2932212) q[3];
sx q[3];
rz(2.3310272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9206813) q[0];
sx q[0];
rz(-2.8784316) q[0];
sx q[0];
rz(2.8780908) q[0];
rz(-1.1031411) q[1];
sx q[1];
rz(-1.3145072) q[1];
sx q[1];
rz(-0.37364328) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6275245) q[0];
sx q[0];
rz(-1.6434945) q[0];
sx q[0];
rz(3.0811653) q[0];
rz(-3.0194026) q[2];
sx q[2];
rz(-0.8475248) q[2];
sx q[2];
rz(-0.10553372) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3195575) q[1];
sx q[1];
rz(-2.4801835) q[1];
sx q[1];
rz(2.103785) q[1];
rz(-pi) q[2];
rz(-2.3832537) q[3];
sx q[3];
rz(-0.74303526) q[3];
sx q[3];
rz(-2.5287927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.11334795) q[2];
sx q[2];
rz(-2.9714606) q[2];
sx q[2];
rz(-2.587758) q[2];
rz(-1.7438186) q[3];
sx q[3];
rz(-2.5388986) q[3];
sx q[3];
rz(0.3127313) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58436191) q[0];
sx q[0];
rz(-1.0065684) q[0];
sx q[0];
rz(1.2297909) q[0];
rz(0.23904414) q[1];
sx q[1];
rz(-1.6300647) q[1];
sx q[1];
rz(-0.30034932) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25324437) q[0];
sx q[0];
rz(-1.7333366) q[0];
sx q[0];
rz(-1.6881315) q[0];
rz(-pi) q[1];
rz(-2.5744252) q[2];
sx q[2];
rz(-0.78710273) q[2];
sx q[2];
rz(1.2445104) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19076599) q[1];
sx q[1];
rz(-3.1263104) q[1];
sx q[1];
rz(-1.5962259) q[1];
rz(-pi) q[2];
rz(3.1170397) q[3];
sx q[3];
rz(-0.85806393) q[3];
sx q[3];
rz(0.36125444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.131669) q[2];
sx q[2];
rz(-1.6841623) q[2];
sx q[2];
rz(2.8966676) q[2];
rz(-2.6217672) q[3];
sx q[3];
rz(-2.2782785) q[3];
sx q[3];
rz(-0.68827099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7898665) q[0];
sx q[0];
rz(-2.7930197) q[0];
sx q[0];
rz(-1.9785731) q[0];
rz(-3.0746958) q[1];
sx q[1];
rz(-1.6480548) q[1];
sx q[1];
rz(2.1287207) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6179919) q[0];
sx q[0];
rz(-0.80192425) q[0];
sx q[0];
rz(-0.84419925) q[0];
rz(-pi) q[1];
rz(0.17872058) q[2];
sx q[2];
rz(-1.0105437) q[2];
sx q[2];
rz(-2.5541039) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0433181) q[1];
sx q[1];
rz(-1.2407082) q[1];
sx q[1];
rz(1.9041474) q[1];
rz(-pi) q[2];
rz(-0.68655218) q[3];
sx q[3];
rz(-1.2408537) q[3];
sx q[3];
rz(0.8984962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.11671994) q[2];
sx q[2];
rz(-0.97815424) q[2];
sx q[2];
rz(-2.8187974) q[2];
rz(-0.60574496) q[3];
sx q[3];
rz(-2.3492458) q[3];
sx q[3];
rz(-2.7927223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054319687) q[0];
sx q[0];
rz(-0.1467341) q[0];
sx q[0];
rz(0.012454575) q[0];
rz(-0.74673486) q[1];
sx q[1];
rz(-0.92339271) q[1];
sx q[1];
rz(2.8616203) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075172193) q[0];
sx q[0];
rz(-0.1233347) q[0];
sx q[0];
rz(-1.146011) q[0];
rz(-pi) q[1];
rz(0.32692636) q[2];
sx q[2];
rz(-1.1571615) q[2];
sx q[2];
rz(1.7978316) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2586883) q[1];
sx q[1];
rz(-1.8778442) q[1];
sx q[1];
rz(0.76489246) q[1];
rz(-pi) q[2];
rz(-0.98484184) q[3];
sx q[3];
rz(-0.42736125) q[3];
sx q[3];
rz(0.40287429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3296457) q[2];
sx q[2];
rz(-0.75355607) q[2];
sx q[2];
rz(2.9296056) q[2];
rz(0.82344615) q[3];
sx q[3];
rz(-1.6971089) q[3];
sx q[3];
rz(-2.8823891) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9987746) q[0];
sx q[0];
rz(-3.0840254) q[0];
sx q[0];
rz(2.4488191) q[0];
rz(-2.5686) q[1];
sx q[1];
rz(-1.3447821) q[1];
sx q[1];
rz(2.7105892) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8739024) q[0];
sx q[0];
rz(-2.5866716) q[0];
sx q[0];
rz(-0.11124994) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5334778) q[2];
sx q[2];
rz(-2.4466142) q[2];
sx q[2];
rz(-2.8335477) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6083518) q[1];
sx q[1];
rz(-2.2679288) q[1];
sx q[1];
rz(-1.9014386) q[1];
rz(-pi) q[2];
rz(2.836976) q[3];
sx q[3];
rz(-1.0416789) q[3];
sx q[3];
rz(-1.2347138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7196322) q[2];
sx q[2];
rz(-0.27406359) q[2];
sx q[2];
rz(-2.5813622) q[2];
rz(-0.46960056) q[3];
sx q[3];
rz(-0.40265366) q[3];
sx q[3];
rz(0.71389055) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2847168) q[0];
sx q[0];
rz(-1.7354043) q[0];
sx q[0];
rz(2.0866557) q[0];
rz(-2.3003385) q[1];
sx q[1];
rz(-1.1019191) q[1];
sx q[1];
rz(-0.046774653) q[1];
rz(1.2011436) q[2];
sx q[2];
rz(-1.2751725) q[2];
sx q[2];
rz(-1.0775492) q[2];
rz(0.026997707) q[3];
sx q[3];
rz(-2.229573) q[3];
sx q[3];
rz(0.98462478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
