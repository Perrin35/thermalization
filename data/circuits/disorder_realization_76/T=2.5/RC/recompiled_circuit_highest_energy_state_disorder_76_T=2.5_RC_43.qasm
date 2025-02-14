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
rz(-0.96867222) q[0];
sx q[0];
rz(2.2157366) q[0];
rz(0.61697382) q[1];
sx q[1];
rz(-0.66520912) q[1];
sx q[1];
rz(-1.8170504) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90841328) q[0];
sx q[0];
rz(-0.97199134) q[0];
sx q[0];
rz(-1.4935059) q[0];
x q[1];
rz(2.5568004) q[2];
sx q[2];
rz(-1.6521786) q[2];
sx q[2];
rz(2.6150102) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.87781964) q[1];
sx q[1];
rz(-0.78697453) q[1];
sx q[1];
rz(-2.469648) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67456986) q[3];
sx q[3];
rz(-0.88817443) q[3];
sx q[3];
rz(-2.5404524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.8806184) q[2];
sx q[2];
rz(-1.3857434) q[2];
sx q[2];
rz(2.3495038) q[2];
rz(-0.875862) q[3];
sx q[3];
rz(-3.0328817) q[3];
sx q[3];
rz(-0.66332269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51139128) q[0];
sx q[0];
rz(-1.317861) q[0];
sx q[0];
rz(-0.9414916) q[0];
rz(-0.9990274) q[1];
sx q[1];
rz(-0.92728725) q[1];
sx q[1];
rz(0.082854465) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0339399) q[0];
sx q[0];
rz(-1.4705188) q[0];
sx q[0];
rz(-2.4058008) q[0];
rz(-pi) q[1];
rz(-0.24225927) q[2];
sx q[2];
rz(-1.9177116) q[2];
sx q[2];
rz(-2.0478947) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8585457) q[1];
sx q[1];
rz(-1.0632005) q[1];
sx q[1];
rz(1.0495484) q[1];
rz(-pi) q[2];
rz(-2.1762455) q[3];
sx q[3];
rz(-2.1248528) q[3];
sx q[3];
rz(-0.29747552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2443709) q[2];
sx q[2];
rz(-1.3542078) q[2];
sx q[2];
rz(1.8841891) q[2];
rz(-0.8463549) q[3];
sx q[3];
rz(-2.9823163) q[3];
sx q[3];
rz(2.4634821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.9427247) q[0];
sx q[0];
rz(-1.4582448) q[0];
sx q[0];
rz(2.9123836) q[0];
rz(2.9275059) q[1];
sx q[1];
rz(-0.87132088) q[1];
sx q[1];
rz(0.3814989) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59807459) q[0];
sx q[0];
rz(-1.2208573) q[0];
sx q[0];
rz(-0.66611163) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3029331) q[2];
sx q[2];
rz(-1.4268954) q[2];
sx q[2];
rz(-3.0001907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2974071) q[1];
sx q[1];
rz(-1.5420037) q[1];
sx q[1];
rz(-1.5178568) q[1];
rz(-pi) q[2];
rz(-2.7945064) q[3];
sx q[3];
rz(-1.575609) q[3];
sx q[3];
rz(-2.6322685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4182338) q[2];
sx q[2];
rz(-0.25622076) q[2];
sx q[2];
rz(0.56824938) q[2];
rz(-1.7226284) q[3];
sx q[3];
rz(-1.7122372) q[3];
sx q[3];
rz(-1.0770146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1133465) q[0];
sx q[0];
rz(-1.1320817) q[0];
sx q[0];
rz(-2.8908253) q[0];
rz(-0.084550683) q[1];
sx q[1];
rz(-2.0471579) q[1];
sx q[1];
rz(-1.9409404) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6687209) q[0];
sx q[0];
rz(-2.3568912) q[0];
sx q[0];
rz(1.0502104) q[0];
rz(-pi) q[1];
rz(2.0617101) q[2];
sx q[2];
rz(-0.49188313) q[2];
sx q[2];
rz(1.5668004) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8507217) q[1];
sx q[1];
rz(-0.71858908) q[1];
sx q[1];
rz(0.1512089) q[1];
x q[2];
rz(-1.0446928) q[3];
sx q[3];
rz(-1.3402437) q[3];
sx q[3];
rz(-1.5073204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67479977) q[2];
sx q[2];
rz(-1.1226706) q[2];
sx q[2];
rz(-1.2910845) q[2];
rz(-0.42168266) q[3];
sx q[3];
rz(-2.6899874) q[3];
sx q[3];
rz(1.565056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6186433) q[0];
sx q[0];
rz(-0.72294253) q[0];
sx q[0];
rz(-1.6408386) q[0];
rz(-0.067642637) q[1];
sx q[1];
rz(-1.5875971) q[1];
sx q[1];
rz(-1.1281475) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85894205) q[0];
sx q[0];
rz(-2.9500486) q[0];
sx q[0];
rz(2.9454977) q[0];
rz(1.0447776) q[2];
sx q[2];
rz(-1.2630208) q[2];
sx q[2];
rz(-0.36164944) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0813898) q[1];
sx q[1];
rz(-0.12559016) q[1];
sx q[1];
rz(1.1274028) q[1];
rz(-pi) q[2];
rz(2.9077072) q[3];
sx q[3];
rz(-1.2986168) q[3];
sx q[3];
rz(0.50258499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.03380123) q[2];
sx q[2];
rz(-1.1148323) q[2];
sx q[2];
rz(-0.6768674) q[2];
rz(0.90773165) q[3];
sx q[3];
rz(-1.2264484) q[3];
sx q[3];
rz(-0.41456732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9229729) q[0];
sx q[0];
rz(-0.75160471) q[0];
sx q[0];
rz(-0.76939097) q[0];
rz(1.2409302) q[1];
sx q[1];
rz(-2.2878094) q[1];
sx q[1];
rz(-0.21533899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.365959) q[0];
sx q[0];
rz(-1.3833519) q[0];
sx q[0];
rz(-2.1523317) q[0];
rz(1.9015997) q[2];
sx q[2];
rz(-1.5546521) q[2];
sx q[2];
rz(0.46260297) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.29500719) q[1];
sx q[1];
rz(-0.93448105) q[1];
sx q[1];
rz(-1.2656487) q[1];
x q[2];
rz(0.39677119) q[3];
sx q[3];
rz(-0.7378398) q[3];
sx q[3];
rz(-2.6529907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.26663366) q[2];
sx q[2];
rz(-1.5077488) q[2];
sx q[2];
rz(-0.61666644) q[2];
rz(-2.941361) q[3];
sx q[3];
rz(-0.73035208) q[3];
sx q[3];
rz(-2.0088137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013833372) q[0];
sx q[0];
rz(-0.56202373) q[0];
sx q[0];
rz(-0.01509893) q[0];
rz(3.0191782) q[1];
sx q[1];
rz(-1.8465123) q[1];
sx q[1];
rz(-1.0282358) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6262344) q[0];
sx q[0];
rz(-1.5152103) q[0];
sx q[0];
rz(-2.4856011) q[0];
rz(-0.9527377) q[2];
sx q[2];
rz(-1.7430804) q[2];
sx q[2];
rz(-0.034644459) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9837357) q[1];
sx q[1];
rz(-2.7089556) q[1];
sx q[1];
rz(-0.067072596) q[1];
rz(-0.31556231) q[3];
sx q[3];
rz(-1.3129819) q[3];
sx q[3];
rz(0.97168789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5369109) q[2];
sx q[2];
rz(-1.2102419) q[2];
sx q[2];
rz(2.6417285) q[2];
rz(0.31575051) q[3];
sx q[3];
rz(-0.67086589) q[3];
sx q[3];
rz(2.1226814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.4278118) q[0];
sx q[0];
rz(-1.918387) q[0];
sx q[0];
rz(2.1977303) q[0];
rz(1.4048514) q[1];
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
rz(-2.0269146) q[0];
sx q[0];
rz(-1.8887547) q[0];
sx q[0];
rz(-1.6284031) q[0];
x q[1];
rz(0.2537639) q[2];
sx q[2];
rz(-1.2813338) q[2];
sx q[2];
rz(2.8004405) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7870175) q[1];
sx q[1];
rz(-1.0598618) q[1];
sx q[1];
rz(-2.4035011) q[1];
rz(-pi) q[2];
rz(1.5783903) q[3];
sx q[3];
rz(-2.6287327) q[3];
sx q[3];
rz(0.22554413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9965808) q[2];
sx q[2];
rz(-1.8700446) q[2];
sx q[2];
rz(-1.5550295) q[2];
rz(2.0874713) q[3];
sx q[3];
rz(-1.8127706) q[3];
sx q[3];
rz(-1.740295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15040511) q[0];
sx q[0];
rz(-1.0973955) q[0];
sx q[0];
rz(-0.64252585) q[0];
rz(-1.7036899) q[1];
sx q[1];
rz(-0.75690401) q[1];
sx q[1];
rz(2.9978851) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7724975) q[0];
sx q[0];
rz(-1.0155627) q[0];
sx q[0];
rz(-1.2124208) q[0];
x q[1];
rz(-0.92387284) q[2];
sx q[2];
rz(-2.0028186) q[2];
sx q[2];
rz(0.96521711) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2255114) q[1];
sx q[1];
rz(-0.71174946) q[1];
sx q[1];
rz(0.41081482) q[1];
rz(2.7034143) q[3];
sx q[3];
rz(-1.2196121) q[3];
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
rz(-2.883834) q[3];
sx q[3];
rz(-0.38193211) q[3];
sx q[3];
rz(-0.8530544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.2951374) q[0];
sx q[0];
rz(-1.4895804) q[0];
sx q[0];
rz(1.9211796) q[0];
rz(-2.659761) q[1];
sx q[1];
rz(-2.316663) q[1];
sx q[1];
rz(-2.2442472) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98678714) q[0];
sx q[0];
rz(-0.53175321) q[0];
sx q[0];
rz(-0.24260862) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8209575) q[2];
sx q[2];
rz(-2.8047987) q[2];
sx q[2];
rz(-1.8048665) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.362295) q[1];
sx q[1];
rz(-0.26302281) q[1];
sx q[1];
rz(-1.3785326) q[1];
x q[2];
rz(-3.0378689) q[3];
sx q[3];
rz(-2.1770766) q[3];
sx q[3];
rz(2.1192354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64409488) q[2];
sx q[2];
rz(-0.92659014) q[2];
sx q[2];
rz(-3.0268055) q[2];
rz(-2.3980906) q[3];
sx q[3];
rz(-1.2657575) q[3];
sx q[3];
rz(-2.6928597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4153618) q[0];
sx q[0];
rz(-0.76114934) q[0];
sx q[0];
rz(2.1834955) q[0];
rz(1.505898) q[1];
sx q[1];
rz(-2.0151357) q[1];
sx q[1];
rz(-1.9912079) q[1];
rz(2.3957797) q[2];
sx q[2];
rz(-2.0796761) q[2];
sx q[2];
rz(2.833119) q[2];
rz(-2.7076247) q[3];
sx q[3];
rz(-1.7229736) q[3];
sx q[3];
rz(2.4255502) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
