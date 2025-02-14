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
rz(0.089601547) q[0];
sx q[0];
rz(2.2360585) q[0];
sx q[0];
rz(9.6146248) q[0];
rz(-0.39659652) q[1];
sx q[1];
rz(-2.8009669) q[1];
sx q[1];
rz(0.61000282) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1413795) q[0];
sx q[0];
rz(-1.3249517) q[0];
sx q[0];
rz(2.8971998) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0135569) q[2];
sx q[2];
rz(-1.4024078) q[2];
sx q[2];
rz(-1.5284644) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.09874092) q[1];
sx q[1];
rz(-2.0804278) q[1];
sx q[1];
rz(2.3364728) q[1];
rz(-pi) q[2];
rz(-1.0759495) q[3];
sx q[3];
rz(-0.54022721) q[3];
sx q[3];
rz(-1.5264508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6923328) q[2];
sx q[2];
rz(-1.924694) q[2];
sx q[2];
rz(-2.8999691) q[2];
rz(0.062945098) q[3];
sx q[3];
rz(-1.9398305) q[3];
sx q[3];
rz(2.3285749) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45074201) q[0];
sx q[0];
rz(-1.8618604) q[0];
sx q[0];
rz(-2.8513841) q[0];
rz(2.8065575) q[1];
sx q[1];
rz(-2.2033043) q[1];
sx q[1];
rz(-0.32523528) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27446628) q[0];
sx q[0];
rz(-0.46733213) q[0];
sx q[0];
rz(0.090666002) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89128691) q[2];
sx q[2];
rz(-1.7789947) q[2];
sx q[2];
rz(3.1058499) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8222386) q[1];
sx q[1];
rz(-1.7387973) q[1];
sx q[1];
rz(0.15138123) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4899204) q[3];
sx q[3];
rz(-2.044994) q[3];
sx q[3];
rz(1.6722566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30514303) q[2];
sx q[2];
rz(-2.2025351) q[2];
sx q[2];
rz(1.1391501) q[2];
rz(0.75974733) q[3];
sx q[3];
rz(-0.097948827) q[3];
sx q[3];
rz(0.37794149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8515795) q[0];
sx q[0];
rz(-0.72999287) q[0];
sx q[0];
rz(-0.83339018) q[0];
rz(0.003412811) q[1];
sx q[1];
rz(-2.623787) q[1];
sx q[1];
rz(1.980967) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7297555) q[0];
sx q[0];
rz(-1.5999536) q[0];
sx q[0];
rz(-1.8957608) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4936475) q[2];
sx q[2];
rz(-1.3908236) q[2];
sx q[2];
rz(2.5729716) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5382446) q[1];
sx q[1];
rz(-1.5711492) q[1];
sx q[1];
rz(-1.3360436) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9594775) q[3];
sx q[3];
rz(-1.210375) q[3];
sx q[3];
rz(0.52626901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87355906) q[2];
sx q[2];
rz(-1.4159091) q[2];
sx q[2];
rz(2.9746919) q[2];
rz(1.9514826) q[3];
sx q[3];
rz(-0.24470617) q[3];
sx q[3];
rz(1.083583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13389182) q[0];
sx q[0];
rz(-1.5465443) q[0];
sx q[0];
rz(-2.290945) q[0];
rz(-0.8029241) q[1];
sx q[1];
rz(-1.9839169) q[1];
sx q[1];
rz(1.1319152) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23861215) q[0];
sx q[0];
rz(-2.065326) q[0];
sx q[0];
rz(2.935845) q[0];
rz(2.1430127) q[2];
sx q[2];
rz(-2.7957186) q[2];
sx q[2];
rz(1.5550176) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.67805659) q[1];
sx q[1];
rz(-0.66036036) q[1];
sx q[1];
rz(-2.6956431) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6904039) q[3];
sx q[3];
rz(-2.3590617) q[3];
sx q[3];
rz(0.53451371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6159281) q[2];
sx q[2];
rz(-1.5027081) q[2];
sx q[2];
rz(-0.28044236) q[2];
rz(-0.40124714) q[3];
sx q[3];
rz(-0.29956996) q[3];
sx q[3];
rz(-1.7846599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2111557) q[0];
sx q[0];
rz(-1.7481952) q[0];
sx q[0];
rz(-2.9537971) q[0];
rz(-0.64741778) q[1];
sx q[1];
rz(-2.1719666) q[1];
sx q[1];
rz(0.59026778) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85344106) q[0];
sx q[0];
rz(-2.7803505) q[0];
sx q[0];
rz(1.4807184) q[0];
x q[1];
rz(-0.36358707) q[2];
sx q[2];
rz(-2.3191904) q[2];
sx q[2];
rz(2.7261312) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5178211) q[1];
sx q[1];
rz(-1.54269) q[1];
sx q[1];
rz(-0.080561056) q[1];
rz(-pi) q[2];
rz(0.26735683) q[3];
sx q[3];
rz(-0.68213576) q[3];
sx q[3];
rz(0.29058829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4831627) q[2];
sx q[2];
rz(-1.1700609) q[2];
sx q[2];
rz(-2.9675193) q[2];
rz(2.6744794) q[3];
sx q[3];
rz(-0.73254782) q[3];
sx q[3];
rz(-2.9113801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9664522) q[0];
sx q[0];
rz(-1.7564961) q[0];
sx q[0];
rz(-1.0873644) q[0];
rz(-0.1611791) q[1];
sx q[1];
rz(-1.9728856) q[1];
sx q[1];
rz(-0.93343121) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95240923) q[0];
sx q[0];
rz(-1.2564474) q[0];
sx q[0];
rz(-1.0745924) q[0];
x q[1];
rz(1.9447787) q[2];
sx q[2];
rz(-1.6829964) q[2];
sx q[2];
rz(-2.2430132) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78466985) q[1];
sx q[1];
rz(-1.9013199) q[1];
sx q[1];
rz(0.4504942) q[1];
rz(2.3087738) q[3];
sx q[3];
rz(-2.3721381) q[3];
sx q[3];
rz(-0.11386816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8072088) q[2];
sx q[2];
rz(-1.7024567) q[2];
sx q[2];
rz(-1.2558233) q[2];
rz(-3.1387591) q[3];
sx q[3];
rz(-0.56003672) q[3];
sx q[3];
rz(-2.475256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46347076) q[0];
sx q[0];
rz(-0.016962873) q[0];
sx q[0];
rz(2.8609138) q[0];
rz(-0.30411389) q[1];
sx q[1];
rz(-2.3801443) q[1];
sx q[1];
rz(-1.4842518) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87697498) q[0];
sx q[0];
rz(-1.5881946) q[0];
sx q[0];
rz(-1.5511284) q[0];
rz(1.7353457) q[2];
sx q[2];
rz(-0.79716792) q[2];
sx q[2];
rz(-2.2556502) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6929984) q[1];
sx q[1];
rz(-1.0585551) q[1];
sx q[1];
rz(2.1650141) q[1];
rz(-pi) q[2];
rz(-1.2400886) q[3];
sx q[3];
rz(-2.7484012) q[3];
sx q[3];
rz(-2.0313944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7228399) q[2];
sx q[2];
rz(-2.3374228) q[2];
sx q[2];
rz(-2.1201102) q[2];
rz(-2.569765) q[3];
sx q[3];
rz(-0.56838667) q[3];
sx q[3];
rz(2.2233326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21917139) q[0];
sx q[0];
rz(-0.89184856) q[0];
sx q[0];
rz(3.0944371) q[0];
rz(-1.7258518) q[1];
sx q[1];
rz(-1.9784617) q[1];
sx q[1];
rz(0.39173752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4120868) q[0];
sx q[0];
rz(-2.0238618) q[0];
sx q[0];
rz(1.6184676) q[0];
rz(-pi) q[1];
rz(-0.18905807) q[2];
sx q[2];
rz(-1.9044276) q[2];
sx q[2];
rz(0.58429532) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89501666) q[1];
sx q[1];
rz(-1.5575174) q[1];
sx q[1];
rz(3.1311656) q[1];
rz(1.3545808) q[3];
sx q[3];
rz(-1.9554227) q[3];
sx q[3];
rz(0.011354488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91531104) q[2];
sx q[2];
rz(-0.87258029) q[2];
sx q[2];
rz(0.20381168) q[2];
rz(-2.5053744) q[3];
sx q[3];
rz(-1.7812984) q[3];
sx q[3];
rz(0.067559592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5485452) q[0];
sx q[0];
rz(-0.021634463) q[0];
sx q[0];
rz(1.1170603) q[0];
rz(1.8695658) q[1];
sx q[1];
rz(-2.2950324) q[1];
sx q[1];
rz(0.55714947) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8654738) q[0];
sx q[0];
rz(-1.1412924) q[0];
sx q[0];
rz(0.19709023) q[0];
rz(2.0821758) q[2];
sx q[2];
rz(-2.9520028) q[2];
sx q[2];
rz(1.0230881) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8073358) q[1];
sx q[1];
rz(-0.42667056) q[1];
sx q[1];
rz(2.2183499) q[1];
rz(0.54160133) q[3];
sx q[3];
rz(-1.3256875) q[3];
sx q[3];
rz(2.4211943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7316526) q[2];
sx q[2];
rz(-0.2356379) q[2];
sx q[2];
rz(1.5059936) q[2];
rz(-2.951494) q[3];
sx q[3];
rz(-0.59571576) q[3];
sx q[3];
rz(0.070040919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6079123) q[0];
sx q[0];
rz(-2.8354225) q[0];
sx q[0];
rz(-0.26468563) q[0];
rz(-2.7429122) q[1];
sx q[1];
rz(-0.52972263) q[1];
sx q[1];
rz(0.20464373) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084572425) q[0];
sx q[0];
rz(-2.1099732) q[0];
sx q[0];
rz(0.23410213) q[0];
rz(-1.9098656) q[2];
sx q[2];
rz(-2.481784) q[2];
sx q[2];
rz(1.8383023) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.8890587) q[1];
sx q[1];
rz(-1.6315347) q[1];
sx q[1];
rz(0.82734682) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7162303) q[3];
sx q[3];
rz(-1.0974636) q[3];
sx q[3];
rz(2.4504599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.22145049) q[2];
sx q[2];
rz(-1.7998989) q[2];
sx q[2];
rz(-1.3915001) q[2];
rz(2.5642388) q[3];
sx q[3];
rz(-0.57830638) q[3];
sx q[3];
rz(-2.3277843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53265424) q[0];
sx q[0];
rz(-1.2464936) q[0];
sx q[0];
rz(-2.0410224) q[0];
rz(-0.34250034) q[1];
sx q[1];
rz(-2.1771912) q[1];
sx q[1];
rz(2.049581) q[1];
rz(-1.0625381) q[2];
sx q[2];
rz(-2.7926805) q[2];
sx q[2];
rz(2.8624848) q[2];
rz(2.2130421) q[3];
sx q[3];
rz(-0.83893574) q[3];
sx q[3];
rz(-0.75225603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
