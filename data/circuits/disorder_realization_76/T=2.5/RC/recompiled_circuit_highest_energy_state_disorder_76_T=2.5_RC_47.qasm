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
rz(0.87585706) q[0];
sx q[0];
rz(2.1729204) q[0];
sx q[0];
rz(7.2090413) q[0];
rz(0.61697382) q[1];
sx q[1];
rz(-0.66520912) q[1];
sx q[1];
rz(1.3245423) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61875859) q[0];
sx q[0];
rz(-1.5069739) q[0];
sx q[0];
rz(0.60019779) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5568004) q[2];
sx q[2];
rz(-1.6521786) q[2];
sx q[2];
rz(-0.52658242) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1086667) q[1];
sx q[1];
rz(-0.98331988) q[1];
sx q[1];
rz(-2.1290178) q[1];
x q[2];
rz(0.76530568) q[3];
sx q[3];
rz(-1.0648784) q[3];
sx q[3];
rz(-1.4368865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2609743) q[2];
sx q[2];
rz(-1.7558492) q[2];
sx q[2];
rz(0.79208881) q[2];
rz(-0.875862) q[3];
sx q[3];
rz(-0.10871092) q[3];
sx q[3];
rz(0.66332269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6302014) q[0];
sx q[0];
rz(-1.317861) q[0];
sx q[0];
rz(-2.200101) q[0];
rz(2.1425653) q[1];
sx q[1];
rz(-0.92728725) q[1];
sx q[1];
rz(0.082854465) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44646406) q[0];
sx q[0];
rz(-0.83954357) q[0];
sx q[0];
rz(1.7056998) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98495667) q[2];
sx q[2];
rz(-0.42030605) q[2];
sx q[2];
rz(2.676385) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99020834) q[1];
sx q[1];
rz(-2.4308204) q[1];
sx q[1];
rz(-0.73020331) q[1];
rz(2.1762455) q[3];
sx q[3];
rz(-1.0167398) q[3];
sx q[3];
rz(-0.29747552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.89722172) q[2];
sx q[2];
rz(-1.3542078) q[2];
sx q[2];
rz(-1.8841891) q[2];
rz(-0.8463549) q[3];
sx q[3];
rz(-0.1592764) q[3];
sx q[3];
rz(0.67811051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-1.9427247) q[0];
sx q[0];
rz(-1.4582448) q[0];
sx q[0];
rz(-2.9123836) q[0];
rz(0.21408679) q[1];
sx q[1];
rz(-0.87132088) q[1];
sx q[1];
rz(-0.3814989) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3840528) q[0];
sx q[0];
rz(-2.4017576) q[0];
sx q[0];
rz(-2.6081041) q[0];
x q[1];
rz(1.3573285) q[2];
sx q[2];
rz(-0.74356438) q[2];
sx q[2];
rz(1.5539757) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.22893045) q[1];
sx q[1];
rz(-0.060256392) q[1];
sx q[1];
rz(1.0723537) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.014147357) q[3];
sx q[3];
rz(-2.7944744) q[3];
sx q[3];
rz(-2.0668168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7233589) q[2];
sx q[2];
rz(-0.25622076) q[2];
sx q[2];
rz(-2.5733433) q[2];
rz(1.7226284) q[3];
sx q[3];
rz(-1.4293554) q[3];
sx q[3];
rz(2.0645781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028246183) q[0];
sx q[0];
rz(-1.1320817) q[0];
sx q[0];
rz(-2.8908253) q[0];
rz(0.084550683) q[1];
sx q[1];
rz(-1.0944347) q[1];
sx q[1];
rz(1.2006522) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4290743) q[0];
sx q[0];
rz(-1.2116644) q[0];
sx q[0];
rz(2.284689) q[0];
x q[1];
rz(1.1293639) q[2];
sx q[2];
rz(-1.3462634) q[2];
sx q[2];
rz(-2.6972636) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6509962) q[1];
sx q[1];
rz(-0.8621434) q[1];
sx q[1];
rz(-1.7017822) q[1];
x q[2];
rz(-2.0968998) q[3];
sx q[3];
rz(-1.3402437) q[3];
sx q[3];
rz(-1.6342722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.67479977) q[2];
sx q[2];
rz(-2.0189221) q[2];
sx q[2];
rz(-1.8505081) q[2];
rz(-2.71991) q[3];
sx q[3];
rz(-0.45160523) q[3];
sx q[3];
rz(-1.5765367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52294937) q[0];
sx q[0];
rz(-0.72294253) q[0];
sx q[0];
rz(-1.6408386) q[0];
rz(3.07395) q[1];
sx q[1];
rz(-1.5875971) q[1];
sx q[1];
rz(2.0134451) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2826506) q[0];
sx q[0];
rz(-0.19154405) q[0];
sx q[0];
rz(0.19609496) q[0];
rz(-pi) q[1];
rz(2.0968151) q[2];
sx q[2];
rz(-1.8785718) q[2];
sx q[2];
rz(-0.36164944) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.060202816) q[1];
sx q[1];
rz(-3.0160025) q[1];
sx q[1];
rz(-1.1274028) q[1];
rz(-1.2913843) q[3];
sx q[3];
rz(-1.7959204) q[3];
sx q[3];
rz(-2.1373419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.03380123) q[2];
sx q[2];
rz(-2.0267603) q[2];
sx q[2];
rz(-0.6768674) q[2];
rz(2.233861) q[3];
sx q[3];
rz(-1.2264484) q[3];
sx q[3];
rz(0.41456732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9229729) q[0];
sx q[0];
rz(-0.75160471) q[0];
sx q[0];
rz(0.76939097) q[0];
rz(-1.2409302) q[1];
sx q[1];
rz(-0.85378328) q[1];
sx q[1];
rz(-0.21533899) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0701987) q[0];
sx q[0];
rz(-0.60766534) q[0];
sx q[0];
rz(-1.2383226) q[0];
rz(-pi) q[1];
rz(-3.1245232) q[2];
sx q[2];
rz(-1.2400377) q[2];
sx q[2];
rz(1.1137373) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0508381) q[1];
sx q[1];
rz(-1.3267446) q[1];
sx q[1];
rz(0.65907101) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2329383) q[3];
sx q[3];
rz(-2.2400899) q[3];
sx q[3];
rz(3.11495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.874959) q[2];
sx q[2];
rz(-1.6338438) q[2];
sx q[2];
rz(-0.61666644) q[2];
rz(-2.941361) q[3];
sx q[3];
rz(-2.4112406) q[3];
sx q[3];
rz(-1.1327789) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1277593) q[0];
sx q[0];
rz(-0.56202373) q[0];
sx q[0];
rz(3.1264937) q[0];
rz(0.12241441) q[1];
sx q[1];
rz(-1.8465123) q[1];
sx q[1];
rz(1.0282358) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1249753) q[0];
sx q[0];
rz(-2.4835971) q[0];
sx q[0];
rz(-3.0506177) q[0];
rz(1.2790658) q[2];
sx q[2];
rz(-0.63858596) q[2];
sx q[2];
rz(-1.8420458) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0575917) q[1];
sx q[1];
rz(-2.0023953) q[1];
sx q[1];
rz(-1.6017385) q[1];
x q[2];
rz(0.31556231) q[3];
sx q[3];
rz(-1.8286108) q[3];
sx q[3];
rz(-2.1699048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5369109) q[2];
sx q[2];
rz(-1.2102419) q[2];
sx q[2];
rz(-2.6417285) q[2];
rz(0.31575051) q[3];
sx q[3];
rz(-2.4707268) q[3];
sx q[3];
rz(-2.1226814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.4278118) q[0];
sx q[0];
rz(-1.2232057) q[0];
sx q[0];
rz(-2.1977303) q[0];
rz(-1.7367412) q[1];
sx q[1];
rz(-1.2662788) q[1];
sx q[1];
rz(2.2199383) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2093344) q[0];
sx q[0];
rz(-2.8186322) q[0];
sx q[0];
rz(-0.17318101) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2537639) q[2];
sx q[2];
rz(-1.8602588) q[2];
sx q[2];
rz(2.8004405) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6347305) q[1];
sx q[1];
rz(-0.94351879) q[1];
sx q[1];
rz(0.92232134) q[1];
rz(-1.5783903) q[3];
sx q[3];
rz(-2.6287327) q[3];
sx q[3];
rz(-0.22554413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1450119) q[2];
sx q[2];
rz(-1.2715481) q[2];
sx q[2];
rz(1.5865631) q[2];
rz(-2.0874713) q[3];
sx q[3];
rz(-1.328822) q[3];
sx q[3];
rz(-1.740295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9911875) q[0];
sx q[0];
rz(-2.0441971) q[0];
sx q[0];
rz(0.64252585) q[0];
rz(-1.7036899) q[1];
sx q[1];
rz(-2.3846886) q[1];
sx q[1];
rz(0.14370758) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7724975) q[0];
sx q[0];
rz(-1.0155627) q[0];
sx q[0];
rz(-1.2124208) q[0];
rz(-pi) q[1];
rz(-2.2238067) q[2];
sx q[2];
rz(-0.76030234) q[2];
sx q[2];
rz(2.0298634) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.66389304) q[1];
sx q[1];
rz(-1.3069005) q[1];
sx q[1];
rz(-2.4724835) q[1];
rz(0.43817839) q[3];
sx q[3];
rz(-1.9219805) q[3];
sx q[3];
rz(-2.8170057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6045427) q[2];
sx q[2];
rz(-1.1979251) q[2];
sx q[2];
rz(2.878888) q[2];
rz(0.25775868) q[3];
sx q[3];
rz(-0.38193211) q[3];
sx q[3];
rz(2.2885382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8464552) q[0];
sx q[0];
rz(-1.4895804) q[0];
sx q[0];
rz(1.9211796) q[0];
rz(0.48183164) q[1];
sx q[1];
rz(-2.316663) q[1];
sx q[1];
rz(-2.2442472) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7677444) q[0];
sx q[0];
rz(-1.4486827) q[0];
sx q[0];
rz(-0.51885817) q[0];
rz(1.4608896) q[2];
sx q[2];
rz(-1.2517923) q[2];
sx q[2];
rz(-1.6750592) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7792977) q[1];
sx q[1];
rz(-0.26302281) q[1];
sx q[1];
rz(1.7630601) q[1];
rz(2.1796024) q[3];
sx q[3];
rz(-1.4856087) q[3];
sx q[3];
rz(-2.5339068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64409488) q[2];
sx q[2];
rz(-2.2150025) q[2];
sx q[2];
rz(-0.11478718) q[2];
rz(0.74350205) q[3];
sx q[3];
rz(-1.8758352) q[3];
sx q[3];
rz(2.6928597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4153618) q[0];
sx q[0];
rz(-2.3804433) q[0];
sx q[0];
rz(-0.95809715) q[0];
rz(1.6356946) q[1];
sx q[1];
rz(-1.1264569) q[1];
sx q[1];
rz(1.1503848) q[1];
rz(-2.4534879) q[2];
sx q[2];
rz(-0.87439121) q[2];
sx q[2];
rz(1.7476358) q[2];
rz(0.34974364) q[3];
sx q[3];
rz(-2.6833224) q[3];
sx q[3];
rz(0.53862017) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
