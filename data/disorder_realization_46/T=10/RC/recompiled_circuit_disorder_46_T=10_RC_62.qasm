OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5503791) q[0];
sx q[0];
rz(-3.1382635) q[0];
sx q[0];
rz(0.34531265) q[0];
rz(1.8058864) q[1];
sx q[1];
rz(3.4808573) q[1];
sx q[1];
rz(9.1453287) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5292408) q[0];
sx q[0];
rz(-1.2623598) q[0];
sx q[0];
rz(-0.2150857) q[0];
rz(1.9136393) q[2];
sx q[2];
rz(-2.3228085) q[2];
sx q[2];
rz(-2.0057099) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.93912032) q[1];
sx q[1];
rz(-1.7608374) q[1];
sx q[1];
rz(3.128016) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0885987) q[3];
sx q[3];
rz(-0.64265673) q[3];
sx q[3];
rz(-1.7794123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44148579) q[2];
sx q[2];
rz(-1.6690648) q[2];
sx q[2];
rz(-2.1842365) q[2];
rz(-0.29933128) q[3];
sx q[3];
rz(-2.744031) q[3];
sx q[3];
rz(0.41199747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9807724) q[0];
sx q[0];
rz(-0.94962025) q[0];
sx q[0];
rz(-0.41369307) q[0];
rz(1.7970239) q[1];
sx q[1];
rz(-2.3584056) q[1];
sx q[1];
rz(2.5059674) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3755075) q[0];
sx q[0];
rz(-1.6151531) q[0];
sx q[0];
rz(-1.7031329) q[0];
rz(-pi) q[1];
rz(0.046819709) q[2];
sx q[2];
rz(-1.4390107) q[2];
sx q[2];
rz(1.5891967) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0483413) q[1];
sx q[1];
rz(-1.8257358) q[1];
sx q[1];
rz(-1.3202207) q[1];
rz(-0.73511519) q[3];
sx q[3];
rz(-2.2921836) q[3];
sx q[3];
rz(-1.848156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3893434) q[2];
sx q[2];
rz(-1.2717335) q[2];
sx q[2];
rz(-0.52865571) q[2];
rz(1.2403437) q[3];
sx q[3];
rz(-0.35959187) q[3];
sx q[3];
rz(-0.24578978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56025958) q[0];
sx q[0];
rz(-1.154705) q[0];
sx q[0];
rz(-2.8833959) q[0];
rz(-1.5064346) q[1];
sx q[1];
rz(-2.5841027) q[1];
sx q[1];
rz(0.43446508) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70179825) q[0];
sx q[0];
rz(-0.5895624) q[0];
sx q[0];
rz(-2.2586285) q[0];
rz(-pi) q[1];
rz(-0.065504727) q[2];
sx q[2];
rz(-2.0993877) q[2];
sx q[2];
rz(-0.62578177) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1669958) q[1];
sx q[1];
rz(-1.3612862) q[1];
sx q[1];
rz(2.6796883) q[1];
x q[2];
rz(-0.86359777) q[3];
sx q[3];
rz(-2.2798385) q[3];
sx q[3];
rz(-1.9883224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7210641) q[2];
sx q[2];
rz(-1.3238182) q[2];
sx q[2];
rz(3.0857962) q[2];
rz(-2.2037286) q[3];
sx q[3];
rz(-0.28356975) q[3];
sx q[3];
rz(-0.24308932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043512251) q[0];
sx q[0];
rz(-2.1976017) q[0];
sx q[0];
rz(3.0010624) q[0];
rz(0.17164104) q[1];
sx q[1];
rz(-1.834603) q[1];
sx q[1];
rz(-0.26352873) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4123654) q[0];
sx q[0];
rz(-1.0050887) q[0];
sx q[0];
rz(1.4501146) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5272335) q[2];
sx q[2];
rz(-0.76459568) q[2];
sx q[2];
rz(-1.489153) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9159689) q[1];
sx q[1];
rz(-1.9150503) q[1];
sx q[1];
rz(2.4179789) q[1];
rz(-pi) q[2];
rz(2.2534091) q[3];
sx q[3];
rz(-0.41999751) q[3];
sx q[3];
rz(2.1995467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0984829) q[2];
sx q[2];
rz(-0.3442328) q[2];
sx q[2];
rz(-1.8919224) q[2];
rz(-1.194687) q[3];
sx q[3];
rz(-2.1134816) q[3];
sx q[3];
rz(0.67888129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93976218) q[0];
sx q[0];
rz(-1.8562466) q[0];
sx q[0];
rz(-0.029065954) q[0];
rz(1.4020231) q[1];
sx q[1];
rz(-0.35750917) q[1];
sx q[1];
rz(0.32863858) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5111361) q[0];
sx q[0];
rz(-1.4689313) q[0];
sx q[0];
rz(1.0826375) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8721696) q[2];
sx q[2];
rz(-1.7413483) q[2];
sx q[2];
rz(1.5829057) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8481816) q[1];
sx q[1];
rz(-1.1383346) q[1];
sx q[1];
rz(-1.6407938) q[1];
rz(-0.7469437) q[3];
sx q[3];
rz(-1.9455823) q[3];
sx q[3];
rz(-2.4657616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0236726) q[2];
sx q[2];
rz(-2.6866044) q[2];
sx q[2];
rz(-2.9809791) q[2];
rz(1.9953856) q[3];
sx q[3];
rz(-1.3137772) q[3];
sx q[3];
rz(2.5207991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49204957) q[0];
sx q[0];
rz(-2.2821125) q[0];
sx q[0];
rz(-0.78654003) q[0];
rz(0.37711626) q[1];
sx q[1];
rz(-2.2943594) q[1];
sx q[1];
rz(0.4424817) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1174283) q[0];
sx q[0];
rz(-0.52485835) q[0];
sx q[0];
rz(-0.79458046) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45583506) q[2];
sx q[2];
rz(-2.5172533) q[2];
sx q[2];
rz(-2.3221743) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6742299) q[1];
sx q[1];
rz(-1.2503887) q[1];
sx q[1];
rz(0.25735374) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2551742) q[3];
sx q[3];
rz(-1.481803) q[3];
sx q[3];
rz(0.86252585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4601712) q[2];
sx q[2];
rz(-2.5632016) q[2];
sx q[2];
rz(2.1703413) q[2];
rz(2.4462637) q[3];
sx q[3];
rz(-0.23614241) q[3];
sx q[3];
rz(-0.42738459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63715315) q[0];
sx q[0];
rz(-1.9313066) q[0];
sx q[0];
rz(2.1321645) q[0];
rz(-2.5908453) q[1];
sx q[1];
rz(-1.2754722) q[1];
sx q[1];
rz(1.1632464) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9775987) q[0];
sx q[0];
rz(-1.7118651) q[0];
sx q[0];
rz(-1.2376357) q[0];
rz(2.8909056) q[2];
sx q[2];
rz(-0.87826585) q[2];
sx q[2];
rz(0.78679774) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35112652) q[1];
sx q[1];
rz(-2.8664221) q[1];
sx q[1];
rz(-0.64063425) q[1];
rz(-0.87485119) q[3];
sx q[3];
rz(-1.3303183) q[3];
sx q[3];
rz(1.2547753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.358868) q[2];
sx q[2];
rz(-1.9903368) q[2];
sx q[2];
rz(-2.3434095) q[2];
rz(2.362137) q[3];
sx q[3];
rz(-0.53648406) q[3];
sx q[3];
rz(-1.9688169) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1786132) q[0];
sx q[0];
rz(-2.6586752) q[0];
sx q[0];
rz(3.0187507) q[0];
rz(0.12610647) q[1];
sx q[1];
rz(-1.6364731) q[1];
sx q[1];
rz(-1.925148) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13623304) q[0];
sx q[0];
rz(-1.1598806) q[0];
sx q[0];
rz(-2.4634325) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4037651) q[2];
sx q[2];
rz(-2.0207496) q[2];
sx q[2];
rz(-3.0917167) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5436732) q[1];
sx q[1];
rz(-1.3733175) q[1];
sx q[1];
rz(-0.39241723) q[1];
rz(-2.8669531) q[3];
sx q[3];
rz(-2.5317149) q[3];
sx q[3];
rz(0.21851893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8806261) q[2];
sx q[2];
rz(-0.61769056) q[2];
sx q[2];
rz(3.0984666) q[2];
rz(-2.96636) q[3];
sx q[3];
rz(-0.8774811) q[3];
sx q[3];
rz(-1.5568679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6373428) q[0];
sx q[0];
rz(-3.0872587) q[0];
sx q[0];
rz(-0.20877008) q[0];
rz(1.4978706) q[1];
sx q[1];
rz(-1.6521963) q[1];
sx q[1];
rz(2.0170905) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1632657) q[0];
sx q[0];
rz(-1.4341337) q[0];
sx q[0];
rz(0.13701963) q[0];
rz(-pi) q[1];
rz(0.71473177) q[2];
sx q[2];
rz(-1.3681612) q[2];
sx q[2];
rz(-1.5243901) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.323656) q[1];
sx q[1];
rz(-1.4374104) q[1];
sx q[1];
rz(-0.0011841983) q[1];
rz(-pi) q[2];
rz(0.48891588) q[3];
sx q[3];
rz(-2.206114) q[3];
sx q[3];
rz(0.92632252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0008529) q[2];
sx q[2];
rz(-2.4177987) q[2];
sx q[2];
rz(0.17803426) q[2];
rz(-1.595165) q[3];
sx q[3];
rz(-1.833257) q[3];
sx q[3];
rz(-0.49062887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35266018) q[0];
sx q[0];
rz(-1.9910318) q[0];
sx q[0];
rz(-2.1066522) q[0];
rz(2.3433698) q[1];
sx q[1];
rz(-2.1612576) q[1];
sx q[1];
rz(-2.9916874) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4193383) q[0];
sx q[0];
rz(-2.1968578) q[0];
sx q[0];
rz(1.1115848) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5871928) q[2];
sx q[2];
rz(-1.8573485) q[2];
sx q[2];
rz(-0.80114844) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.80430154) q[1];
sx q[1];
rz(-1.1608487) q[1];
sx q[1];
rz(1.0211584) q[1];
rz(0.30431872) q[3];
sx q[3];
rz(-0.70239866) q[3];
sx q[3];
rz(-1.9866634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4118816) q[2];
sx q[2];
rz(-1.231266) q[2];
sx q[2];
rz(0.40851545) q[2];
rz(0.92489964) q[3];
sx q[3];
rz(-2.3291406) q[3];
sx q[3];
rz(-2.6598721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9020486) q[0];
sx q[0];
rz(-1.5938546) q[0];
sx q[0];
rz(-1.6123733) q[0];
rz(-1.7655903) q[1];
sx q[1];
rz(-1.9648432) q[1];
sx q[1];
rz(1.2479938) q[1];
rz(-1.8953163) q[2];
sx q[2];
rz(-0.80249716) q[2];
sx q[2];
rz(-1.3182166) q[2];
rz(-0.33384791) q[3];
sx q[3];
rz(-0.6469938) q[3];
sx q[3];
rz(-2.7713431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];