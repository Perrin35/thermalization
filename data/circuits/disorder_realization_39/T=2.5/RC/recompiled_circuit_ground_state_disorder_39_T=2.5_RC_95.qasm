OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7118536) q[0];
sx q[0];
rz(-1.1450333) q[0];
sx q[0];
rz(1.0166919) q[0];
rz(-0.84438762) q[1];
sx q[1];
rz(-0.8165741) q[1];
sx q[1];
rz(-2.5752697) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1745464) q[0];
sx q[0];
rz(-0.61515795) q[0];
sx q[0];
rz(1.7394203) q[0];
rz(-pi) q[1];
rz(-0.077458642) q[2];
sx q[2];
rz(-1.3734387) q[2];
sx q[2];
rz(0.058930581) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6366292) q[1];
sx q[1];
rz(-2.7111021) q[1];
sx q[1];
rz(-2.9477547) q[1];
rz(-1.7773966) q[3];
sx q[3];
rz(-0.99541621) q[3];
sx q[3];
rz(0.30295886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6301959) q[2];
sx q[2];
rz(-2.2679195) q[2];
sx q[2];
rz(-2.013618) q[2];
rz(1.5026211) q[3];
sx q[3];
rz(-0.84533826) q[3];
sx q[3];
rz(-0.42475167) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2750435) q[0];
sx q[0];
rz(-0.4568704) q[0];
sx q[0];
rz(-3.1006815) q[0];
rz(0.54967898) q[1];
sx q[1];
rz(-0.37893852) q[1];
sx q[1];
rz(-0.080370195) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59118862) q[0];
sx q[0];
rz(-2.4192717) q[0];
sx q[0];
rz(3.1352311) q[0];
rz(-0.54687107) q[2];
sx q[2];
rz(-2.139922) q[2];
sx q[2];
rz(2.6324181) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6192542) q[1];
sx q[1];
rz(-1.0864746) q[1];
sx q[1];
rz(-0.062141499) q[1];
x q[2];
rz(-1.9527797) q[3];
sx q[3];
rz(-0.63900164) q[3];
sx q[3];
rz(-1.5572214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4139159) q[2];
sx q[2];
rz(-2.2985986) q[2];
sx q[2];
rz(-0.28398871) q[2];
rz(-1.5208288) q[3];
sx q[3];
rz(-1.5512543) q[3];
sx q[3];
rz(-2.3782597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82336998) q[0];
sx q[0];
rz(-2.9871873) q[0];
sx q[0];
rz(0.41686091) q[0];
rz(2.2082224) q[1];
sx q[1];
rz(-0.88637543) q[1];
sx q[1];
rz(-2.0534168) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25167187) q[0];
sx q[0];
rz(-1.3265189) q[0];
sx q[0];
rz(2.1273462) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6112664) q[2];
sx q[2];
rz(-1.010561) q[2];
sx q[2];
rz(-2.3793085) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8510906) q[1];
sx q[1];
rz(-0.41176418) q[1];
sx q[1];
rz(0.74300933) q[1];
x q[2];
rz(-0.21622323) q[3];
sx q[3];
rz(-2.8475886) q[3];
sx q[3];
rz(0.8651498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0205959) q[2];
sx q[2];
rz(-0.72046295) q[2];
sx q[2];
rz(0.32568112) q[2];
rz(2.7459512) q[3];
sx q[3];
rz(-0.49042693) q[3];
sx q[3];
rz(0.71780786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2242103) q[0];
sx q[0];
rz(-2.2069187) q[0];
sx q[0];
rz(0.14064661) q[0];
rz(-0.39528254) q[1];
sx q[1];
rz(-1.544516) q[1];
sx q[1];
rz(0.43221727) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.062586322) q[0];
sx q[0];
rz(-0.64029988) q[0];
sx q[0];
rz(2.2878134) q[0];
rz(-pi) q[1];
rz(-1.725196) q[2];
sx q[2];
rz(-2.6302166) q[2];
sx q[2];
rz(1.3647773) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2905457) q[1];
sx q[1];
rz(-1.6509027) q[1];
sx q[1];
rz(-0.27388957) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3234576) q[3];
sx q[3];
rz(-0.77518089) q[3];
sx q[3];
rz(-1.1213746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.846659) q[2];
sx q[2];
rz(-2.0294971) q[2];
sx q[2];
rz(-0.46889949) q[2];
rz(0.15639671) q[3];
sx q[3];
rz(-2.1726435) q[3];
sx q[3];
rz(1.8739353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073931996) q[0];
sx q[0];
rz(-1.1336552) q[0];
sx q[0];
rz(0.78039783) q[0];
rz(2.228915) q[1];
sx q[1];
rz(-0.78881216) q[1];
sx q[1];
rz(-1.7524293) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71043832) q[0];
sx q[0];
rz(-1.9214993) q[0];
sx q[0];
rz(-1.2855269) q[0];
rz(-0.19540968) q[2];
sx q[2];
rz(-0.59048803) q[2];
sx q[2];
rz(2.5974642) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9174172) q[1];
sx q[1];
rz(-0.97748102) q[1];
sx q[1];
rz(1.7043056) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0205998) q[3];
sx q[3];
rz(-0.76348272) q[3];
sx q[3];
rz(-0.17624763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3461561) q[2];
sx q[2];
rz(-0.6554335) q[2];
sx q[2];
rz(-0.1795086) q[2];
rz(1.4795715) q[3];
sx q[3];
rz(-0.69412762) q[3];
sx q[3];
rz(-0.0050541335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(0.96106225) q[0];
sx q[0];
rz(-1.5164277) q[0];
sx q[0];
rz(-1.6777212) q[0];
rz(-0.013710984) q[1];
sx q[1];
rz(-1.4669712) q[1];
sx q[1];
rz(0.94508583) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4714519) q[0];
sx q[0];
rz(-1.2218804) q[0];
sx q[0];
rz(2.3405911) q[0];
rz(-pi) q[1];
rz(2.5588972) q[2];
sx q[2];
rz(-2.7740363) q[2];
sx q[2];
rz(2.3876397) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2234563) q[1];
sx q[1];
rz(-2.0386749) q[1];
sx q[1];
rz(2.0418057) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1626445) q[3];
sx q[3];
rz(-1.1148387) q[3];
sx q[3];
rz(0.78612529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1549687) q[2];
sx q[2];
rz(-2.0782317) q[2];
sx q[2];
rz(-1.0490136) q[2];
rz(1.5341885) q[3];
sx q[3];
rz(-2.1395855) q[3];
sx q[3];
rz(1.775942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55117115) q[0];
sx q[0];
rz(-0.48749247) q[0];
sx q[0];
rz(-2.3792939) q[0];
rz(-0.47362348) q[1];
sx q[1];
rz(-1.381424) q[1];
sx q[1];
rz(0.90829888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.565641) q[0];
sx q[0];
rz(-1.5759908) q[0];
sx q[0];
rz(-2.3495102) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90877675) q[2];
sx q[2];
rz(-0.27784608) q[2];
sx q[2];
rz(1.2706333) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1181284) q[1];
sx q[1];
rz(-2.9042517) q[1];
sx q[1];
rz(-1.508839) q[1];
x q[2];
rz(-1.9929816) q[3];
sx q[3];
rz(-0.030638846) q[3];
sx q[3];
rz(-2.2922761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0198387) q[2];
sx q[2];
rz(-2.8826931) q[2];
sx q[2];
rz(0.961595) q[2];
rz(-2.614295) q[3];
sx q[3];
rz(-1.7617825) q[3];
sx q[3];
rz(0.56078792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90802646) q[0];
sx q[0];
rz(-2.5416424) q[0];
sx q[0];
rz(0.34616923) q[0];
rz(-1.2973805) q[1];
sx q[1];
rz(-2.4751414) q[1];
sx q[1];
rz(2.5328439) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93701836) q[0];
sx q[0];
rz(-1.4868344) q[0];
sx q[0];
rz(2.2473826) q[0];
rz(-pi) q[1];
rz(2.8738451) q[2];
sx q[2];
rz(-1.7742947) q[2];
sx q[2];
rz(-0.45058077) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.94683719) q[1];
sx q[1];
rz(-1.5018592) q[1];
sx q[1];
rz(-0.26534715) q[1];
x q[2];
rz(0.99842803) q[3];
sx q[3];
rz(-2.2791036) q[3];
sx q[3];
rz(0.64841333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79494563) q[2];
sx q[2];
rz(-2.5865159) q[2];
sx q[2];
rz(0.65183276) q[2];
rz(0.73976222) q[3];
sx q[3];
rz(-2.0756105) q[3];
sx q[3];
rz(-2.2447926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.5480492) q[0];
sx q[0];
rz(-2.6683922) q[0];
sx q[0];
rz(-0.59762534) q[0];
rz(-1.301544) q[1];
sx q[1];
rz(-1.5664682) q[1];
sx q[1];
rz(-2.8155933) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38446389) q[0];
sx q[0];
rz(-1.281503) q[0];
sx q[0];
rz(2.3209653) q[0];
rz(-pi) q[1];
rz(-1.3131485) q[2];
sx q[2];
rz(-1.8882053) q[2];
sx q[2];
rz(2.7820415) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.9218775) q[1];
sx q[1];
rz(-1.3120717) q[1];
sx q[1];
rz(2.2398276) q[1];
rz(-1.4526618) q[3];
sx q[3];
rz(-1.6836327) q[3];
sx q[3];
rz(-2.0510587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43805435) q[2];
sx q[2];
rz(-1.9376829) q[2];
sx q[2];
rz(-0.34029141) q[2];
rz(-1.5002286) q[3];
sx q[3];
rz(-0.54268018) q[3];
sx q[3];
rz(2.5435508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.930645) q[0];
sx q[0];
rz(-1.1177381) q[0];
sx q[0];
rz(-3.0979544) q[0];
rz(2.4234407) q[1];
sx q[1];
rz(-1.3190045) q[1];
sx q[1];
rz(1.7125548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0008416) q[0];
sx q[0];
rz(-1.4006097) q[0];
sx q[0];
rz(-0.19077459) q[0];
rz(2.9731644) q[2];
sx q[2];
rz(-2.456924) q[2];
sx q[2];
rz(-1.1845961) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.90718118) q[1];
sx q[1];
rz(-2.6446807) q[1];
sx q[1];
rz(0.18053825) q[1];
rz(-1.1943071) q[3];
sx q[3];
rz(-1.912279) q[3];
sx q[3];
rz(-1.3545413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4622197) q[2];
sx q[2];
rz(-0.93928176) q[2];
sx q[2];
rz(-0.77793724) q[2];
rz(-2.098162) q[3];
sx q[3];
rz(-1.5305887) q[3];
sx q[3];
rz(1.9061609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0935681) q[0];
sx q[0];
rz(-1.0014191) q[0];
sx q[0];
rz(0.77145664) q[0];
rz(1.363516) q[1];
sx q[1];
rz(-1.9693146) q[1];
sx q[1];
rz(1.6065425) q[1];
rz(2.5566348) q[2];
sx q[2];
rz(-0.81333209) q[2];
sx q[2];
rz(1.3438136) q[2];
rz(-0.89911581) q[3];
sx q[3];
rz(-1.9439144) q[3];
sx q[3];
rz(1.6685566) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
