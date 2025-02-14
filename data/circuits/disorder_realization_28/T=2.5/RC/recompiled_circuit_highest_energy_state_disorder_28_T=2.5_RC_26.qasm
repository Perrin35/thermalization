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
rz(0.75594354) q[0];
sx q[0];
rz(3.4721952) q[0];
sx q[0];
rz(9.1520206) q[0];
rz(2.7859712) q[1];
sx q[1];
rz(-2.1115117) q[1];
sx q[1];
rz(-1.570809) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.49424) q[0];
sx q[0];
rz(-1.6675416) q[0];
sx q[0];
rz(0.072150196) q[0];
rz(-1.2414957) q[2];
sx q[2];
rz(-1.1277756) q[2];
sx q[2];
rz(1.7274513) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.65793234) q[1];
sx q[1];
rz(-2.2495553) q[1];
sx q[1];
rz(1.5163877) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4727888) q[3];
sx q[3];
rz(-1.9379172) q[3];
sx q[3];
rz(-2.844732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.97751272) q[2];
sx q[2];
rz(-1.7089607) q[2];
sx q[2];
rz(2.254503) q[2];
rz(-1.2648434) q[3];
sx q[3];
rz(-1.0586459) q[3];
sx q[3];
rz(3.1177055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7510659) q[0];
sx q[0];
rz(-0.85705119) q[0];
sx q[0];
rz(2.8476025) q[0];
rz(-1.3514306) q[1];
sx q[1];
rz(-2.0452812) q[1];
sx q[1];
rz(-1.9757804) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6751926) q[0];
sx q[0];
rz(-1.8972016) q[0];
sx q[0];
rz(-0.2194039) q[0];
rz(1.4977156) q[2];
sx q[2];
rz(-1.5433307) q[2];
sx q[2];
rz(2.0171201) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.15951751) q[1];
sx q[1];
rz(-2.0802092) q[1];
sx q[1];
rz(-1.0647091) q[1];
x q[2];
rz(0.56880672) q[3];
sx q[3];
rz(-1.8173462) q[3];
sx q[3];
rz(-2.5914471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.68842781) q[2];
sx q[2];
rz(-0.63899779) q[2];
sx q[2];
rz(1.2028018) q[2];
rz(1.5996251) q[3];
sx q[3];
rz(-1.8034233) q[3];
sx q[3];
rz(2.8560824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0315392) q[0];
sx q[0];
rz(-3.0634614) q[0];
sx q[0];
rz(1.9670638) q[0];
rz(0.9749167) q[1];
sx q[1];
rz(-2.2623623) q[1];
sx q[1];
rz(2.6851795) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5147091) q[0];
sx q[0];
rz(-2.868342) q[0];
sx q[0];
rz(1.68863) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24791294) q[2];
sx q[2];
rz(-1.1968024) q[2];
sx q[2];
rz(-3.0097112) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.2172981) q[1];
sx q[1];
rz(-0.84554377) q[1];
sx q[1];
rz(0.78332232) q[1];
x q[2];
rz(2.1022845) q[3];
sx q[3];
rz(-2.3520497) q[3];
sx q[3];
rz(-1.2293775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5259033) q[2];
sx q[2];
rz(-2.758785) q[2];
sx q[2];
rz(2.3846386) q[2];
rz(-3.1225539) q[3];
sx q[3];
rz(-1.2913387) q[3];
sx q[3];
rz(2.3058057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2794936) q[0];
sx q[0];
rz(-2.4308496) q[0];
sx q[0];
rz(2.1918462) q[0];
rz(2.9277335) q[1];
sx q[1];
rz(-2.2016134) q[1];
sx q[1];
rz(-2.5423539) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99357641) q[0];
sx q[0];
rz(-0.92088078) q[0];
sx q[0];
rz(0.17493576) q[0];
x q[1];
rz(-1.2798115) q[2];
sx q[2];
rz(-2.6070234) q[2];
sx q[2];
rz(-0.9415516) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.039700198) q[1];
sx q[1];
rz(-1.3624316) q[1];
sx q[1];
rz(1.031967) q[1];
rz(-pi) q[2];
rz(0.20885373) q[3];
sx q[3];
rz(-0.59126544) q[3];
sx q[3];
rz(2.8610736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5102101) q[2];
sx q[2];
rz(-1.5082794) q[2];
sx q[2];
rz(0.82376662) q[2];
rz(1.0907762) q[3];
sx q[3];
rz(-2.3407276) q[3];
sx q[3];
rz(-2.4763079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3983696) q[0];
sx q[0];
rz(-2.7581765) q[0];
sx q[0];
rz(-2.1552591) q[0];
rz(0.24364722) q[1];
sx q[1];
rz(-1.8758834) q[1];
sx q[1];
rz(-2.9820014) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.900478) q[0];
sx q[0];
rz(-2.3474446) q[0];
sx q[0];
rz(-2.0179124) q[0];
rz(-2.1417732) q[2];
sx q[2];
rz(-2.2717064) q[2];
sx q[2];
rz(2.4510405) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3978093) q[1];
sx q[1];
rz(-1.127217) q[1];
sx q[1];
rz(-3.0438782) q[1];
rz(-pi) q[2];
rz(0.87017228) q[3];
sx q[3];
rz(-1.3425832) q[3];
sx q[3];
rz(0.96503497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.39449447) q[2];
sx q[2];
rz(-1.2622086) q[2];
sx q[2];
rz(-1.1355642) q[2];
rz(-0.46843946) q[3];
sx q[3];
rz(-1.9197542) q[3];
sx q[3];
rz(0.92456094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0202476) q[0];
sx q[0];
rz(-2.0609042) q[0];
sx q[0];
rz(-0.24170804) q[0];
rz(1.3555591) q[1];
sx q[1];
rz(-0.87239289) q[1];
sx q[1];
rz(-2.2579069) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51118219) q[0];
sx q[0];
rz(-0.97994655) q[0];
sx q[0];
rz(2.1686337) q[0];
rz(-pi) q[1];
rz(2.8293351) q[2];
sx q[2];
rz(-1.3259282) q[2];
sx q[2];
rz(-1.6409525) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7979398) q[1];
sx q[1];
rz(-0.48544297) q[1];
sx q[1];
rz(1.9992827) q[1];
rz(-3.1413636) q[3];
sx q[3];
rz(-1.9817686) q[3];
sx q[3];
rz(-2.4966893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.22150618) q[2];
sx q[2];
rz(-0.93044996) q[2];
sx q[2];
rz(2.7152756) q[2];
rz(0.91655556) q[3];
sx q[3];
rz(-1.5896268) q[3];
sx q[3];
rz(-1.1030654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78974462) q[0];
sx q[0];
rz(-1.9639356) q[0];
sx q[0];
rz(1.1258997) q[0];
rz(-3.0806091) q[1];
sx q[1];
rz(-2.1911502) q[1];
sx q[1];
rz(-2.1263863) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.489129) q[0];
sx q[0];
rz(-1.154107) q[0];
sx q[0];
rz(0.66834895) q[0];
rz(-pi) q[1];
rz(0.073748579) q[2];
sx q[2];
rz(-1.908386) q[2];
sx q[2];
rz(0.13843564) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8422994) q[1];
sx q[1];
rz(-1.1720578) q[1];
sx q[1];
rz(0.066940424) q[1];
x q[2];
rz(-1.6105846) q[3];
sx q[3];
rz(-1.6841751) q[3];
sx q[3];
rz(-1.6802406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8288377) q[2];
sx q[2];
rz(-1.2042896) q[2];
sx q[2];
rz(0.22605669) q[2];
rz(2.1894646) q[3];
sx q[3];
rz(-3.003037) q[3];
sx q[3];
rz(0.8086732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011618135) q[0];
sx q[0];
rz(-1.3148774) q[0];
sx q[0];
rz(2.8011978) q[0];
rz(0.099523425) q[1];
sx q[1];
rz(-2.0390022) q[1];
sx q[1];
rz(1.5460825) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9700546) q[0];
sx q[0];
rz(-2.4405789) q[0];
sx q[0];
rz(-2.6215977) q[0];
rz(-3.0853188) q[2];
sx q[2];
rz(-1.2985126) q[2];
sx q[2];
rz(0.19993609) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9118285) q[1];
sx q[1];
rz(-0.10783261) q[1];
sx q[1];
rz(-0.39984881) q[1];
rz(1.4813878) q[3];
sx q[3];
rz(-2.0242888) q[3];
sx q[3];
rz(0.88624533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.16326363) q[2];
sx q[2];
rz(-1.949387) q[2];
sx q[2];
rz(-1.0015798) q[2];
rz(1.352365) q[3];
sx q[3];
rz(-0.47836256) q[3];
sx q[3];
rz(1.9549687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1773862) q[0];
sx q[0];
rz(-1.0564251) q[0];
sx q[0];
rz(0.00038432234) q[0];
rz(-2.7035825) q[1];
sx q[1];
rz(-1.507746) q[1];
sx q[1];
rz(0.45723525) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23822902) q[0];
sx q[0];
rz(-1.5212802) q[0];
sx q[0];
rz(-1.7286506) q[0];
x q[1];
rz(0.57827823) q[2];
sx q[2];
rz(-2.1077029) q[2];
sx q[2];
rz(-2.0202877) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.25343382) q[1];
sx q[1];
rz(-2.7256084) q[1];
sx q[1];
rz(2.5640954) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5874279) q[3];
sx q[3];
rz(-0.57798741) q[3];
sx q[3];
rz(-0.91855861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.75927258) q[2];
sx q[2];
rz(-2.1542532) q[2];
sx q[2];
rz(-1.6299204) q[2];
rz(3.0509389) q[3];
sx q[3];
rz(-2.1001308) q[3];
sx q[3];
rz(-0.67971984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7882305) q[0];
sx q[0];
rz(-1.3557949) q[0];
sx q[0];
rz(-2.8613388) q[0];
rz(-1.3837815) q[1];
sx q[1];
rz(-0.46418142) q[1];
sx q[1];
rz(-1.2253449) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89577755) q[0];
sx q[0];
rz(-1.7208463) q[0];
sx q[0];
rz(1.4729985) q[0];
rz(-1.046243) q[2];
sx q[2];
rz(-0.75511564) q[2];
sx q[2];
rz(-0.81449997) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8333421) q[1];
sx q[1];
rz(-2.3099515) q[1];
sx q[1];
rz(-1.9019446) q[1];
x q[2];
rz(-1.2810542) q[3];
sx q[3];
rz(-1.8612766) q[3];
sx q[3];
rz(1.8673459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.19259351) q[2];
sx q[2];
rz(-1.5263564) q[2];
sx q[2];
rz(0.68142778) q[2];
rz(2.0415908) q[3];
sx q[3];
rz(-0.93183485) q[3];
sx q[3];
rz(-1.3070235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.6901907) q[0];
sx q[0];
rz(-2.207353) q[0];
sx q[0];
rz(1.9718715) q[0];
rz(-2.7611217) q[1];
sx q[1];
rz(-2.4741551) q[1];
sx q[1];
rz(-2.9173775) q[1];
rz(2.2718471) q[2];
sx q[2];
rz(-1.7746584) q[2];
sx q[2];
rz(-0.45523582) q[2];
rz(-2.7093665) q[3];
sx q[3];
rz(-2.1189011) q[3];
sx q[3];
rz(-0.44174474) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
