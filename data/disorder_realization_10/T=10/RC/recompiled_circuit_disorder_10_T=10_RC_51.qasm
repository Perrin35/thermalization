OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(-1.9987885) q[0];
sx q[0];
rz(-1.9300652) q[0];
rz(6.05655) q[1];
sx q[1];
rz(1.5645138) q[1];
sx q[1];
rz(5.984879) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5814712) q[0];
sx q[0];
rz(-1.5157962) q[0];
sx q[0];
rz(-0.66530692) q[0];
rz(-pi) q[1];
rz(3.0478893) q[2];
sx q[2];
rz(-1.9022577) q[2];
sx q[2];
rz(-1.1790438) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6023941) q[1];
sx q[1];
rz(-1.0326003) q[1];
sx q[1];
rz(0.56166517) q[1];
rz(-pi) q[2];
rz(1.467642) q[3];
sx q[3];
rz(-1.6543596) q[3];
sx q[3];
rz(0.32675693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6478708) q[2];
sx q[2];
rz(-1.2080668) q[2];
sx q[2];
rz(-2.1477264) q[2];
rz(0.99938756) q[3];
sx q[3];
rz(-1.9013654) q[3];
sx q[3];
rz(-0.58888155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64269972) q[0];
sx q[0];
rz(-1.9357341) q[0];
sx q[0];
rz(3.1233741) q[0];
rz(2.3253564) q[1];
sx q[1];
rz(-2.1111592) q[1];
sx q[1];
rz(2.6699064) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1175849) q[0];
sx q[0];
rz(-1.4837259) q[0];
sx q[0];
rz(-2.9136806) q[0];
rz(-pi) q[1];
rz(2.6623146) q[2];
sx q[2];
rz(-0.55715484) q[2];
sx q[2];
rz(0.011475871) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3018803) q[1];
sx q[1];
rz(-1.0186968) q[1];
sx q[1];
rz(-1.2549972) q[1];
rz(-pi) q[2];
rz(2.9845893) q[3];
sx q[3];
rz(-0.96884851) q[3];
sx q[3];
rz(0.70037819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54962426) q[2];
sx q[2];
rz(-1.9886118) q[2];
sx q[2];
rz(-2.1726051) q[2];
rz(-2.5668868) q[3];
sx q[3];
rz(-2.59022) q[3];
sx q[3];
rz(-1.0415174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7085003) q[0];
sx q[0];
rz(-2.0860465) q[0];
sx q[0];
rz(-0.18181268) q[0];
rz(1.1026985) q[1];
sx q[1];
rz(-1.5010553) q[1];
sx q[1];
rz(-1.4556494) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91711125) q[0];
sx q[0];
rz(-0.79229504) q[0];
sx q[0];
rz(-2.2729421) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4007912) q[2];
sx q[2];
rz(-1.2503137) q[2];
sx q[2];
rz(-1.4893116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.599217) q[1];
sx q[1];
rz(-1.5812751) q[1];
sx q[1];
rz(-0.12565617) q[1];
rz(0.62526838) q[3];
sx q[3];
rz(-1.6571952) q[3];
sx q[3];
rz(0.82739917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2640947) q[2];
sx q[2];
rz(-1.654518) q[2];
sx q[2];
rz(-2.1739615) q[2];
rz(2.4140221) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(-2.9038866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2847292) q[0];
sx q[0];
rz(-0.52607042) q[0];
sx q[0];
rz(-0.63823429) q[0];
rz(-1.1278641) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(1.2329873) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15105948) q[0];
sx q[0];
rz(-2.6418243) q[0];
sx q[0];
rz(-2.997056) q[0];
x q[1];
rz(-0.74027503) q[2];
sx q[2];
rz(-2.4057655) q[2];
sx q[2];
rz(-2.1528113) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2265046) q[1];
sx q[1];
rz(-2.9595032) q[1];
sx q[1];
rz(-1.7712797) q[1];
x q[2];
rz(0.87807699) q[3];
sx q[3];
rz(-1.3948166) q[3];
sx q[3];
rz(0.39719492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91810742) q[2];
sx q[2];
rz(-1.2780259) q[2];
sx q[2];
rz(0.81400648) q[2];
rz(-2.0984086) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(-1.1842747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84213132) q[0];
sx q[0];
rz(-1.786754) q[0];
sx q[0];
rz(0.89170757) q[0];
rz(1.2437598) q[1];
sx q[1];
rz(-1.7638821) q[1];
sx q[1];
rz(-2.9290501) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5355797) q[0];
sx q[0];
rz(-1.5870464) q[0];
sx q[0];
rz(-1.5910801) q[0];
rz(3.1031102) q[2];
sx q[2];
rz(-2.1750692) q[2];
sx q[2];
rz(-2.2272) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7674539) q[1];
sx q[1];
rz(-2.1530188) q[1];
sx q[1];
rz(-2.2351082) q[1];
rz(2.0810633) q[3];
sx q[3];
rz(-0.79499309) q[3];
sx q[3];
rz(1.0802964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1084958) q[2];
sx q[2];
rz(-2.1704845) q[2];
sx q[2];
rz(2.6203716) q[2];
rz(1.7565578) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(-1.8732171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8591156) q[0];
sx q[0];
rz(-2.2127667) q[0];
sx q[0];
rz(2.916472) q[0];
rz(1.7865932) q[1];
sx q[1];
rz(-2.1332108) q[1];
sx q[1];
rz(-2.7640142) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21340428) q[0];
sx q[0];
rz(-1.8367935) q[0];
sx q[0];
rz(-3.1187952) q[0];
x q[1];
rz(-0.16212459) q[2];
sx q[2];
rz(-0.29106489) q[2];
sx q[2];
rz(-1.7578917) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.40124711) q[1];
sx q[1];
rz(-2.0013072) q[1];
sx q[1];
rz(0.11534782) q[1];
rz(-pi) q[2];
rz(1.6013006) q[3];
sx q[3];
rz(-0.77611938) q[3];
sx q[3];
rz(-3.0200849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1569415) q[2];
sx q[2];
rz(-0.7876544) q[2];
sx q[2];
rz(0.48103508) q[2];
rz(-2.7379819) q[3];
sx q[3];
rz(-2.1026881) q[3];
sx q[3];
rz(-0.31479442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.0890546) q[0];
sx q[0];
rz(-2.5367694) q[0];
sx q[0];
rz(-2.9470434) q[0];
rz(0.21952195) q[1];
sx q[1];
rz(-1.6794645) q[1];
sx q[1];
rz(2.887168) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1628111) q[0];
sx q[0];
rz(-0.33919507) q[0];
sx q[0];
rz(1.2374872) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9585602) q[2];
sx q[2];
rz(-1.4441635) q[2];
sx q[2];
rz(2.3938092) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2563045) q[1];
sx q[1];
rz(-1.2997775) q[1];
sx q[1];
rz(-0.46896743) q[1];
x q[2];
rz(1.1021348) q[3];
sx q[3];
rz(-1.9797167) q[3];
sx q[3];
rz(-1.8917781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1640132) q[2];
sx q[2];
rz(-1.7501202) q[2];
sx q[2];
rz(1.8257726) q[2];
rz(-0.87604648) q[3];
sx q[3];
rz(-3.0026569) q[3];
sx q[3];
rz(2.1379437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2106237) q[0];
sx q[0];
rz(-0.3759149) q[0];
sx q[0];
rz(1.4550495) q[0];
rz(-2.3176106) q[1];
sx q[1];
rz(-2.1897557) q[1];
sx q[1];
rz(1.5751858) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9109089) q[0];
sx q[0];
rz(-0.77373234) q[0];
sx q[0];
rz(-0.45549972) q[0];
rz(-1.8640395) q[2];
sx q[2];
rz(-0.55985057) q[2];
sx q[2];
rz(-2.6960889) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7745061) q[1];
sx q[1];
rz(-1.204406) q[1];
sx q[1];
rz(2.7584502) q[1];
x q[2];
rz(0.50623399) q[3];
sx q[3];
rz(-1.1971548) q[3];
sx q[3];
rz(0.95844275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87551293) q[2];
sx q[2];
rz(-1.2337039) q[2];
sx q[2];
rz(2.8640462) q[2];
rz(-1.4510441) q[3];
sx q[3];
rz(-2.6896559) q[3];
sx q[3];
rz(-2.1267166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16185109) q[0];
sx q[0];
rz(-2.6888872) q[0];
sx q[0];
rz(1.6850527) q[0];
rz(-2.5121571) q[1];
sx q[1];
rz(-1.1673735) q[1];
sx q[1];
rz(1.1368407) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6900078) q[0];
sx q[0];
rz(-1.0613872) q[0];
sx q[0];
rz(2.8103229) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37461899) q[2];
sx q[2];
rz(-1.1748474) q[2];
sx q[2];
rz(1.3499201) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6362308) q[1];
sx q[1];
rz(-0.33126918) q[1];
sx q[1];
rz(1.2760217) q[1];
x q[2];
rz(1.9433446) q[3];
sx q[3];
rz(-2.2714943) q[3];
sx q[3];
rz(2.8230599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4920766) q[2];
sx q[2];
rz(-1.7121544) q[2];
sx q[2];
rz(-2.1006404) q[2];
rz(0.0020290931) q[3];
sx q[3];
rz(-2.7414331) q[3];
sx q[3];
rz(0.31203976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8390389) q[0];
sx q[0];
rz(-0.22452393) q[0];
sx q[0];
rz(-2.1955406) q[0];
rz(-0.91167766) q[1];
sx q[1];
rz(-1.9263575) q[1];
sx q[1];
rz(2.5295703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48869041) q[0];
sx q[0];
rz(-2.2757747) q[0];
sx q[0];
rz(-0.61625723) q[0];
rz(-pi) q[1];
rz(0.34531784) q[2];
sx q[2];
rz(-1.3441663) q[2];
sx q[2];
rz(2.741284) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1336085) q[1];
sx q[1];
rz(-2.4719279) q[1];
sx q[1];
rz(2.942251) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6066783) q[3];
sx q[3];
rz(-1.8992918) q[3];
sx q[3];
rz(1.2319777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0845906) q[2];
sx q[2];
rz(-0.64208639) q[2];
sx q[2];
rz(-2.4882312) q[2];
rz(2.7907794) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(-0.70070926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6012797) q[0];
sx q[0];
rz(-1.606034) q[0];
sx q[0];
rz(0.10869797) q[0];
rz(-0.75469771) q[1];
sx q[1];
rz(-1.3194059) q[1];
sx q[1];
rz(-1.5059765) q[1];
rz(-2.4738612) q[2];
sx q[2];
rz(-0.79955352) q[2];
sx q[2];
rz(-0.59895589) q[2];
rz(-0.084372088) q[3];
sx q[3];
rz(-1.8934688) q[3];
sx q[3];
rz(-2.7432751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
