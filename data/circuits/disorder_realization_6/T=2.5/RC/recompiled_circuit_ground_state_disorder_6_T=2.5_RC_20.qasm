OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4226469) q[0];
sx q[0];
rz(-0.5889686) q[0];
sx q[0];
rz(0.11864057) q[0];
rz(-0.9912107) q[1];
sx q[1];
rz(-2.0999496) q[1];
sx q[1];
rz(-0.51377327) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.936718) q[0];
sx q[0];
rz(-0.74130171) q[0];
sx q[0];
rz(-0.78420774) q[0];
rz(-pi) q[1];
rz(-1.5085601) q[2];
sx q[2];
rz(-2.0449491) q[2];
sx q[2];
rz(1.7218423) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.61600191) q[1];
sx q[1];
rz(-0.39037986) q[1];
sx q[1];
rz(-2.5486773) q[1];
rz(-pi) q[2];
rz(-0.051187201) q[3];
sx q[3];
rz(-2.5489106) q[3];
sx q[3];
rz(-0.41646233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.674268) q[2];
sx q[2];
rz(-1.6215723) q[2];
sx q[2];
rz(2.8094214) q[2];
rz(0.25003555) q[3];
sx q[3];
rz(-1.9146405) q[3];
sx q[3];
rz(1.8488319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2501204) q[0];
sx q[0];
rz(-1.253506) q[0];
sx q[0];
rz(2.6943595) q[0];
rz(-1.0294634) q[1];
sx q[1];
rz(-2.5599458) q[1];
sx q[1];
rz(-3.0404125) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6529236) q[0];
sx q[0];
rz(-2.5258668) q[0];
sx q[0];
rz(2.4020345) q[0];
rz(-2.9415628) q[2];
sx q[2];
rz(-2.1406271) q[2];
sx q[2];
rz(2.9404158) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.236233) q[1];
sx q[1];
rz(-1.836685) q[1];
sx q[1];
rz(-2.7120525) q[1];
x q[2];
rz(0.11240837) q[3];
sx q[3];
rz(-0.93494697) q[3];
sx q[3];
rz(-0.67316717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0589361) q[2];
sx q[2];
rz(-0.75581789) q[2];
sx q[2];
rz(-1.5400881) q[2];
rz(-0.61795175) q[3];
sx q[3];
rz(-0.58303419) q[3];
sx q[3];
rz(1.1624973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59442941) q[0];
sx q[0];
rz(-1.0048486) q[0];
sx q[0];
rz(2.6255703) q[0];
rz(-2.0891321) q[1];
sx q[1];
rz(-0.92603374) q[1];
sx q[1];
rz(-1.9693536) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2502278) q[0];
sx q[0];
rz(-2.7173882) q[0];
sx q[0];
rz(-3.1147309) q[0];
x q[1];
rz(-1.6555696) q[2];
sx q[2];
rz(-0.71226487) q[2];
sx q[2];
rz(2.149947) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0111085) q[1];
sx q[1];
rz(-1.7569033) q[1];
sx q[1];
rz(-2.0687201) q[1];
rz(-0.80023685) q[3];
sx q[3];
rz(-1.5808148) q[3];
sx q[3];
rz(-3.0992233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.27511328) q[2];
sx q[2];
rz(-1.2629513) q[2];
sx q[2];
rz(-2.4647554) q[2];
rz(0.46946851) q[3];
sx q[3];
rz(-0.89478409) q[3];
sx q[3];
rz(1.9278056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6770099) q[0];
sx q[0];
rz(-1.0053758) q[0];
sx q[0];
rz(1.9418035) q[0];
rz(-2.5489573) q[1];
sx q[1];
rz(-2.0694144) q[1];
sx q[1];
rz(1.4592272) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30577393) q[0];
sx q[0];
rz(-2.523382) q[0];
sx q[0];
rz(-0.34247663) q[0];
rz(-0.24942579) q[2];
sx q[2];
rz(-0.81455961) q[2];
sx q[2];
rz(0.047175353) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6886518) q[1];
sx q[1];
rz(-2.4600852) q[1];
sx q[1];
rz(1.9863434) q[1];
rz(-pi) q[2];
rz(-1.191889) q[3];
sx q[3];
rz(-0.92273308) q[3];
sx q[3];
rz(-2.5398382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6478641) q[2];
sx q[2];
rz(-1.0961327) q[2];
sx q[2];
rz(2.4596821) q[2];
rz(0.41334263) q[3];
sx q[3];
rz(-1.6091434) q[3];
sx q[3];
rz(1.1558862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2979564) q[0];
sx q[0];
rz(-1.6652668) q[0];
sx q[0];
rz(0.27808878) q[0];
rz(-2.2372712) q[1];
sx q[1];
rz(-1.4537289) q[1];
sx q[1];
rz(2.1542737) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.322418) q[0];
sx q[0];
rz(-1.8083982) q[0];
sx q[0];
rz(2.1850719) q[0];
rz(-pi) q[1];
rz(-0.1369517) q[2];
sx q[2];
rz(-1.5429351) q[2];
sx q[2];
rz(-0.34649039) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5161572) q[1];
sx q[1];
rz(-1.181681) q[1];
sx q[1];
rz(-0.066246943) q[1];
rz(2.0102536) q[3];
sx q[3];
rz(-2.5601697) q[3];
sx q[3];
rz(1.9245337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1272993) q[2];
sx q[2];
rz(-0.96724808) q[2];
sx q[2];
rz(2.3822752) q[2];
rz(-1.6589818) q[3];
sx q[3];
rz(-1.9448152) q[3];
sx q[3];
rz(0.35879859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9229014) q[0];
sx q[0];
rz(-0.8194812) q[0];
sx q[0];
rz(0.6915834) q[0];
rz(2.2870731) q[1];
sx q[1];
rz(-0.71047345) q[1];
sx q[1];
rz(-2.7202594) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74752676) q[0];
sx q[0];
rz(-2.0025064) q[0];
sx q[0];
rz(0.075448087) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26701228) q[2];
sx q[2];
rz(-1.8521554) q[2];
sx q[2];
rz(-2.6002778) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.2273352) q[1];
sx q[1];
rz(-0.37082878) q[1];
sx q[1];
rz(-0.23359681) q[1];
x q[2];
rz(1.297703) q[3];
sx q[3];
rz(-0.98326245) q[3];
sx q[3];
rz(1.3466101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2104346) q[2];
sx q[2];
rz(-1.6633818) q[2];
sx q[2];
rz(2.9388536) q[2];
rz(1.5984795) q[3];
sx q[3];
rz(-2.8212382) q[3];
sx q[3];
rz(-1.6725484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8572674) q[0];
sx q[0];
rz(-1.4284644) q[0];
sx q[0];
rz(2.7427234) q[0];
rz(1.162989) q[1];
sx q[1];
rz(-2.3154924) q[1];
sx q[1];
rz(0.2805447) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6847072) q[0];
sx q[0];
rz(-2.5734757) q[0];
sx q[0];
rz(-0.76596188) q[0];
rz(-pi) q[1];
x q[1];
rz(2.248204) q[2];
sx q[2];
rz(-1.9147493) q[2];
sx q[2];
rz(-2.266573) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.043634) q[1];
sx q[1];
rz(-2.6065738) q[1];
sx q[1];
rz(-2.19997) q[1];
x q[2];
rz(-0.97784247) q[3];
sx q[3];
rz(-0.8675608) q[3];
sx q[3];
rz(0.11444005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.82105381) q[2];
sx q[2];
rz(-2.1619449) q[2];
sx q[2];
rz(1.4873571) q[2];
rz(2.1982543) q[3];
sx q[3];
rz(-0.92883674) q[3];
sx q[3];
rz(1.9802035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5281552) q[0];
sx q[0];
rz(-1.4035839) q[0];
sx q[0];
rz(-0.50233895) q[0];
rz(2.595064) q[1];
sx q[1];
rz(-1.2556475) q[1];
sx q[1];
rz(-2.6796403) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1621409) q[0];
sx q[0];
rz(-0.58425036) q[0];
sx q[0];
rz(-1.0229848) q[0];
rz(-pi) q[1];
rz(-1.5080323) q[2];
sx q[2];
rz(-2.4893005) q[2];
sx q[2];
rz(1.4029274) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4224902) q[1];
sx q[1];
rz(-1.6797069) q[1];
sx q[1];
rz(1.9816573) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5835207) q[3];
sx q[3];
rz(-2.1012348) q[3];
sx q[3];
rz(-0.96543559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.168557) q[2];
sx q[2];
rz(-0.23739561) q[2];
sx q[2];
rz(0.17042223) q[2];
rz(0.45525822) q[3];
sx q[3];
rz(-1.5452496) q[3];
sx q[3];
rz(-2.4236603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1667204) q[0];
sx q[0];
rz(-1.5019187) q[0];
sx q[0];
rz(0.18044743) q[0];
rz(-2.6249053) q[1];
sx q[1];
rz(-1.5645212) q[1];
sx q[1];
rz(2.7076941) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93709263) q[0];
sx q[0];
rz(-0.70815101) q[0];
sx q[0];
rz(2.2592179) q[0];
rz(-pi) q[1];
rz(0.34738587) q[2];
sx q[2];
rz(-1.1466031) q[2];
sx q[2];
rz(3.1372084) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.89131195) q[1];
sx q[1];
rz(-1.4509333) q[1];
sx q[1];
rz(-0.93617546) q[1];
x q[2];
rz(1.9187886) q[3];
sx q[3];
rz(-2.4204614) q[3];
sx q[3];
rz(0.43681555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.58069289) q[2];
sx q[2];
rz(-0.43792024) q[2];
sx q[2];
rz(-1.8221347) q[2];
rz(-2.8699919) q[3];
sx q[3];
rz(-1.3180132) q[3];
sx q[3];
rz(1.1118838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44529799) q[0];
sx q[0];
rz(-0.91067186) q[0];
sx q[0];
rz(-1.4060422) q[0];
rz(1.7156853) q[1];
sx q[1];
rz(-2.0366663) q[1];
sx q[1];
rz(1.8941194) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6794327) q[0];
sx q[0];
rz(-1.832167) q[0];
sx q[0];
rz(1.0699468) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.054385238) q[2];
sx q[2];
rz(-1.453389) q[2];
sx q[2];
rz(-0.94737999) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0669603) q[1];
sx q[1];
rz(-1.6116643) q[1];
sx q[1];
rz(-0.34217477) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91478552) q[3];
sx q[3];
rz(-2.088123) q[3];
sx q[3];
rz(-3.0539258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4282816) q[2];
sx q[2];
rz(-2.9094628) q[2];
sx q[2];
rz(-3.0408119) q[2];
rz(-3.1320324) q[3];
sx q[3];
rz(-1.5567501) q[3];
sx q[3];
rz(-1.5338219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99209256) q[0];
sx q[0];
rz(-1.3652353) q[0];
sx q[0];
rz(1.9944763) q[0];
rz(0.59303444) q[1];
sx q[1];
rz(-1.7094163) q[1];
sx q[1];
rz(-0.97547668) q[1];
rz(-0.070746919) q[2];
sx q[2];
rz(-1.1732709) q[2];
sx q[2];
rz(2.7931961) q[2];
rz(1.5531874) q[3];
sx q[3];
rz(-1.7918158) q[3];
sx q[3];
rz(0.50504897) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
