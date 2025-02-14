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
rz(2.7409878) q[0];
sx q[0];
rz(2.7317943) q[0];
sx q[0];
rz(8.3538342) q[0];
rz(0.61869705) q[1];
sx q[1];
rz(-2.5735811) q[1];
sx q[1];
rz(-0.84363371) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0742357) q[0];
sx q[0];
rz(-1.6003185) q[0];
sx q[0];
rz(1.3062472) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8461692) q[2];
sx q[2];
rz(-1.43338) q[2];
sx q[2];
rz(-0.58770056) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2045793) q[1];
sx q[1];
rz(-0.76256231) q[1];
sx q[1];
rz(2.049905) q[1];
rz(0.47445221) q[3];
sx q[3];
rz(-2.0182924) q[3];
sx q[3];
rz(-1.3407269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7817276) q[2];
sx q[2];
rz(-0.83911506) q[2];
sx q[2];
rz(-0.016999379) q[2];
rz(1.4012236) q[3];
sx q[3];
rz(-0.56178105) q[3];
sx q[3];
rz(-1.9248272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.568999) q[0];
sx q[0];
rz(-1.8091135) q[0];
sx q[0];
rz(-2.3542985) q[0];
rz(0.17838082) q[1];
sx q[1];
rz(-1.9352103) q[1];
sx q[1];
rz(-1.23752) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4263692) q[0];
sx q[0];
rz(-2.5551717) q[0];
sx q[0];
rz(-0.1118025) q[0];
rz(-pi) q[1];
rz(-0.95713116) q[2];
sx q[2];
rz(-1.8120013) q[2];
sx q[2];
rz(-1.7231736) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5893277) q[1];
sx q[1];
rz(-2.5922212) q[1];
sx q[1];
rz(-2.2087847) q[1];
x q[2];
rz(-1.848189) q[3];
sx q[3];
rz(-1.9964868) q[3];
sx q[3];
rz(2.7077146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.42147192) q[2];
sx q[2];
rz(-1.4482435) q[2];
sx q[2];
rz(-0.063701542) q[2];
rz(1.2740159) q[3];
sx q[3];
rz(-2.2985022) q[3];
sx q[3];
rz(1.564285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0434697) q[0];
sx q[0];
rz(-1.948057) q[0];
sx q[0];
rz(0.7435588) q[0];
rz(-0.45383635) q[1];
sx q[1];
rz(-1.791714) q[1];
sx q[1];
rz(-0.41890621) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1322238) q[0];
sx q[0];
rz(-2.0133349) q[0];
sx q[0];
rz(2.924463) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30861993) q[2];
sx q[2];
rz(-0.86008027) q[2];
sx q[2];
rz(-0.98352226) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.031098628) q[1];
sx q[1];
rz(-0.95606177) q[1];
sx q[1];
rz(-0.35009274) q[1];
x q[2];
rz(-2.2499372) q[3];
sx q[3];
rz(-1.4640995) q[3];
sx q[3];
rz(-1.5859491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0539187) q[2];
sx q[2];
rz(-2.3232465) q[2];
sx q[2];
rz(1.9107001) q[2];
rz(0.538921) q[3];
sx q[3];
rz(-1.475622) q[3];
sx q[3];
rz(1.6787136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26381275) q[0];
sx q[0];
rz(-0.75436622) q[0];
sx q[0];
rz(0.60212773) q[0];
rz(-1.4471794) q[1];
sx q[1];
rz(-0.68232957) q[1];
sx q[1];
rz(2.2785861) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9842152) q[0];
sx q[0];
rz(-2.4113305) q[0];
sx q[0];
rz(-2.9709816) q[0];
rz(0.30089897) q[2];
sx q[2];
rz(-1.931012) q[2];
sx q[2];
rz(0.27733251) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80490404) q[1];
sx q[1];
rz(-0.64283481) q[1];
sx q[1];
rz(1.1240608) q[1];
x q[2];
rz(-1.9163777) q[3];
sx q[3];
rz(-2.3732819) q[3];
sx q[3];
rz(2.5490724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0756695) q[2];
sx q[2];
rz(-1.1752335) q[2];
sx q[2];
rz(-0.90265957) q[2];
rz(-0.81501189) q[3];
sx q[3];
rz(-0.60324001) q[3];
sx q[3];
rz(-2.7284315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0937423) q[0];
sx q[0];
rz(-0.047690064) q[0];
sx q[0];
rz(-0.93511859) q[0];
rz(-1.6085666) q[1];
sx q[1];
rz(-2.8159499) q[1];
sx q[1];
rz(-2.8757222) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1810721) q[0];
sx q[0];
rz(-1.4645394) q[0];
sx q[0];
rz(-1.8828859) q[0];
rz(-1.8498672) q[2];
sx q[2];
rz(-0.84013591) q[2];
sx q[2];
rz(1.0026635) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7794687) q[1];
sx q[1];
rz(-1.9196132) q[1];
sx q[1];
rz(-1.7415943) q[1];
rz(-0.10281296) q[3];
sx q[3];
rz(-1.8342557) q[3];
sx q[3];
rz(0.78373945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.461146) q[2];
sx q[2];
rz(-1.5004044) q[2];
sx q[2];
rz(-0.45073304) q[2];
rz(1.9510673) q[3];
sx q[3];
rz(-0.7014941) q[3];
sx q[3];
rz(-0.43959555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9615237) q[0];
sx q[0];
rz(-2.838205) q[0];
sx q[0];
rz(3.0971089) q[0];
rz(-2.8728409) q[1];
sx q[1];
rz(-0.93695295) q[1];
sx q[1];
rz(-2.7106947) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92545729) q[0];
sx q[0];
rz(-1.8325984) q[0];
sx q[0];
rz(2.5039704) q[0];
rz(-2.1304279) q[2];
sx q[2];
rz(-0.75087386) q[2];
sx q[2];
rz(-0.73926413) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0468586) q[1];
sx q[1];
rz(-2.2206535) q[1];
sx q[1];
rz(-2.6296494) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66029064) q[3];
sx q[3];
rz(-1.5662282) q[3];
sx q[3];
rz(-2.3167603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1162794) q[2];
sx q[2];
rz(-2.7677324) q[2];
sx q[2];
rz(2.5171793) q[2];
rz(2.0404909) q[3];
sx q[3];
rz(-0.81239429) q[3];
sx q[3];
rz(2.1541434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2295912) q[0];
sx q[0];
rz(-2.0822552) q[0];
sx q[0];
rz(2.4578995) q[0];
rz(1.4423485) q[1];
sx q[1];
rz(-2.3577299) q[1];
sx q[1];
rz(-0.20194617) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1037052) q[0];
sx q[0];
rz(-1.5695086) q[0];
sx q[0];
rz(-1.8989858) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3181386) q[2];
sx q[2];
rz(-1.5353893) q[2];
sx q[2];
rz(3.0454783) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6579288) q[1];
sx q[1];
rz(-2.6158608) q[1];
sx q[1];
rz(-0.13327285) q[1];
x q[2];
rz(-2.6838164) q[3];
sx q[3];
rz(-1.5321863) q[3];
sx q[3];
rz(-0.82070551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.086229) q[2];
sx q[2];
rz(-1.7217041) q[2];
sx q[2];
rz(-1.8602547) q[2];
rz(2.4751439) q[3];
sx q[3];
rz(-2.0446348) q[3];
sx q[3];
rz(0.54653978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1619038) q[0];
sx q[0];
rz(-0.55210102) q[0];
sx q[0];
rz(-2.6276278) q[0];
rz(-2.1886096) q[1];
sx q[1];
rz(-0.62478462) q[1];
sx q[1];
rz(-1.1187339) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3866769) q[0];
sx q[0];
rz(-2.710452) q[0];
sx q[0];
rz(-2.9453631) q[0];
rz(-0.4301901) q[2];
sx q[2];
rz(-0.83900827) q[2];
sx q[2];
rz(-2.6553287) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2988904) q[1];
sx q[1];
rz(-1.7849079) q[1];
sx q[1];
rz(-1.9886809) q[1];
x q[2];
rz(-1.9268613) q[3];
sx q[3];
rz(-2.0907932) q[3];
sx q[3];
rz(2.4442152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.755456) q[2];
sx q[2];
rz(-2.7894207) q[2];
sx q[2];
rz(0.90292162) q[2];
rz(1.6821945) q[3];
sx q[3];
rz(-1.3420339) q[3];
sx q[3];
rz(-0.82227796) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4626386) q[0];
sx q[0];
rz(-2.8215388) q[0];
sx q[0];
rz(-1.1902887) q[0];
rz(2.8207488) q[1];
sx q[1];
rz(-1.4535934) q[1];
sx q[1];
rz(1.920059) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10968835) q[0];
sx q[0];
rz(-1.3656989) q[0];
sx q[0];
rz(1.431365) q[0];
rz(-pi) q[1];
rz(1.2668816) q[2];
sx q[2];
rz(-0.7579782) q[2];
sx q[2];
rz(2.9044819) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1561191) q[1];
sx q[1];
rz(-1.3669984) q[1];
sx q[1];
rz(-0.76038313) q[1];
rz(1.753958) q[3];
sx q[3];
rz(-1.7622928) q[3];
sx q[3];
rz(-0.93246704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90423501) q[2];
sx q[2];
rz(-0.15494896) q[2];
sx q[2];
rz(2.8083189) q[2];
rz(1.0171558) q[3];
sx q[3];
rz(-2.1867496) q[3];
sx q[3];
rz(0.3199544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79621133) q[0];
sx q[0];
rz(-0.51207241) q[0];
sx q[0];
rz(0.91878015) q[0];
rz(-0.19110075) q[1];
sx q[1];
rz(-0.42593503) q[1];
sx q[1];
rz(2.3474615) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0988172) q[0];
sx q[0];
rz(-2.2198027) q[0];
sx q[0];
rz(-2.2521583) q[0];
rz(-pi) q[1];
rz(-2.7808519) q[2];
sx q[2];
rz(-1.3731925) q[2];
sx q[2];
rz(0.2195356) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0550784) q[1];
sx q[1];
rz(-0.76755133) q[1];
sx q[1];
rz(-1.1805736) q[1];
x q[2];
rz(0.32148949) q[3];
sx q[3];
rz(-1.3703385) q[3];
sx q[3];
rz(-0.050762477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.11067757) q[2];
sx q[2];
rz(-2.5258749) q[2];
sx q[2];
rz(-2.3827629) q[2];
rz(-2.2714254) q[3];
sx q[3];
rz(-0.78799677) q[3];
sx q[3];
rz(-1.0677392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1398685) q[0];
sx q[0];
rz(-1.4586466) q[0];
sx q[0];
rz(1.9139105) q[0];
rz(0.43791804) q[1];
sx q[1];
rz(-1.4358078) q[1];
sx q[1];
rz(-1.5461071) q[1];
rz(-0.6952473) q[2];
sx q[2];
rz(-0.60374478) q[2];
sx q[2];
rz(0.20460621) q[2];
rz(-2.6043456) q[3];
sx q[3];
rz(-0.8575079) q[3];
sx q[3];
rz(-1.433123) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
