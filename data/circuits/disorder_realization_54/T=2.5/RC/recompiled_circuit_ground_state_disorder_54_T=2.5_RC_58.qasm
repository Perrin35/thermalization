OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0889283) q[0];
sx q[0];
rz(5.1533617) q[0];
sx q[0];
rz(12.552153) q[0];
rz(0.05834236) q[1];
sx q[1];
rz(-2.3925233) q[1];
sx q[1];
rz(-2.5383389) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8612807) q[0];
sx q[0];
rz(-2.2076108) q[0];
sx q[0];
rz(2.0122819) q[0];
rz(-pi) q[1];
rz(3.0809359) q[2];
sx q[2];
rz(-0.91691518) q[2];
sx q[2];
rz(-0.72055179) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3989714) q[1];
sx q[1];
rz(-1.4150029) q[1];
sx q[1];
rz(-2.4473518) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7885277) q[3];
sx q[3];
rz(-0.7184295) q[3];
sx q[3];
rz(-1.4361385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.76838905) q[2];
sx q[2];
rz(-1.6037805) q[2];
sx q[2];
rz(1.0962037) q[2];
rz(1.0783892) q[3];
sx q[3];
rz(-0.53467852) q[3];
sx q[3];
rz(2.6064742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32644367) q[0];
sx q[0];
rz(-1.0468227) q[0];
sx q[0];
rz(-0.4775508) q[0];
rz(-2.1304456) q[1];
sx q[1];
rz(-2.5944581) q[1];
sx q[1];
rz(-2.0571041) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80637425) q[0];
sx q[0];
rz(-1.9382297) q[0];
sx q[0];
rz(3.0159877) q[0];
x q[1];
rz(1.162503) q[2];
sx q[2];
rz(-1.5442241) q[2];
sx q[2];
rz(2.0808737) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0140927) q[1];
sx q[1];
rz(-1.0261826) q[1];
sx q[1];
rz(1.7782393) q[1];
rz(0.22900692) q[3];
sx q[3];
rz(-2.3888616) q[3];
sx q[3];
rz(-1.5040042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.61887211) q[2];
sx q[2];
rz(-2.7174157) q[2];
sx q[2];
rz(1.0901394) q[2];
rz(2.5777396) q[3];
sx q[3];
rz(-2.4520051) q[3];
sx q[3];
rz(-0.032698154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3801333) q[0];
sx q[0];
rz(-0.021012336) q[0];
sx q[0];
rz(1.6047961) q[0];
rz(1.8236632) q[1];
sx q[1];
rz(-2.047796) q[1];
sx q[1];
rz(-0.39753786) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56016541) q[0];
sx q[0];
rz(-0.6869964) q[0];
sx q[0];
rz(1.9182179) q[0];
rz(0.2961646) q[2];
sx q[2];
rz(-1.0131179) q[2];
sx q[2];
rz(-0.22685834) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2469625) q[1];
sx q[1];
rz(-1.3161462) q[1];
sx q[1];
rz(2.2130284) q[1];
rz(-pi) q[2];
rz(0.52513772) q[3];
sx q[3];
rz(-1.1236089) q[3];
sx q[3];
rz(2.1839328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90618187) q[2];
sx q[2];
rz(-1.1778888) q[2];
sx q[2];
rz(0.57514352) q[2];
rz(-2.4681674) q[3];
sx q[3];
rz(-1.6005102) q[3];
sx q[3];
rz(-1.5967691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49066082) q[0];
sx q[0];
rz(-2.0199825) q[0];
sx q[0];
rz(-0.90890539) q[0];
rz(-0.017223651) q[1];
sx q[1];
rz(-0.52809087) q[1];
sx q[1];
rz(0.012103279) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7229268) q[0];
sx q[0];
rz(-1.8932492) q[0];
sx q[0];
rz(2.5315383) q[0];
rz(-0.66476314) q[2];
sx q[2];
rz(-2.3741907) q[2];
sx q[2];
rz(0.27945159) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7722655) q[1];
sx q[1];
rz(-0.64913926) q[1];
sx q[1];
rz(1.9499798) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5254667) q[3];
sx q[3];
rz(-1.6974276) q[3];
sx q[3];
rz(1.3172729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.072731344) q[2];
sx q[2];
rz(-2.8673745) q[2];
sx q[2];
rz(0.38632986) q[2];
rz(0.38935152) q[3];
sx q[3];
rz(-1.5334305) q[3];
sx q[3];
rz(-2.1336335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7858081) q[0];
sx q[0];
rz(-2.4171827) q[0];
sx q[0];
rz(-2.8177148) q[0];
rz(0.38662275) q[1];
sx q[1];
rz(-1.7522248) q[1];
sx q[1];
rz(1.5547543) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6232658) q[0];
sx q[0];
rz(-2.4072106) q[0];
sx q[0];
rz(3.0452888) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98734314) q[2];
sx q[2];
rz(-0.61827055) q[2];
sx q[2];
rz(1.5702919) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.388098) q[1];
sx q[1];
rz(-0.83067229) q[1];
sx q[1];
rz(-2.4945033) q[1];
x q[2];
rz(2.6324248) q[3];
sx q[3];
rz(-2.2751791) q[3];
sx q[3];
rz(2.1455163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3891478) q[2];
sx q[2];
rz(-2.2989595) q[2];
sx q[2];
rz(-0.064099224) q[2];
rz(0.60424232) q[3];
sx q[3];
rz(-1.7490381) q[3];
sx q[3];
rz(-1.909168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9391249) q[0];
sx q[0];
rz(-0.61102837) q[0];
sx q[0];
rz(-0.87919277) q[0];
rz(0.35036707) q[1];
sx q[1];
rz(-1.8354974) q[1];
sx q[1];
rz(2.9894357) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6317196) q[0];
sx q[0];
rz(-1.841396) q[0];
sx q[0];
rz(2.4218049) q[0];
rz(0.19292508) q[2];
sx q[2];
rz(-2.0461651) q[2];
sx q[2];
rz(1.0032636) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2895246) q[1];
sx q[1];
rz(-1.3686603) q[1];
sx q[1];
rz(2.7006671) q[1];
rz(0.71546666) q[3];
sx q[3];
rz(-2.5757534) q[3];
sx q[3];
rz(2.6239606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0396314) q[2];
sx q[2];
rz(-2.8122718) q[2];
sx q[2];
rz(0.24001089) q[2];
rz(-0.65252423) q[3];
sx q[3];
rz(-2.0034761) q[3];
sx q[3];
rz(1.2751381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0535102) q[0];
sx q[0];
rz(-1.9048012) q[0];
sx q[0];
rz(-2.2834593) q[0];
rz(1.3872604) q[1];
sx q[1];
rz(-0.45448449) q[1];
sx q[1];
rz(-0.33285704) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7031606) q[0];
sx q[0];
rz(-1.0045719) q[0];
sx q[0];
rz(0.37545183) q[0];
rz(-pi) q[1];
rz(0.22631876) q[2];
sx q[2];
rz(-2.8853325) q[2];
sx q[2];
rz(1.2761436) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7499125) q[1];
sx q[1];
rz(-1.9894162) q[1];
sx q[1];
rz(1.3454354) q[1];
x q[2];
rz(1.9234571) q[3];
sx q[3];
rz(-0.94649345) q[3];
sx q[3];
rz(-1.2144292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3420928) q[2];
sx q[2];
rz(-2.9436593) q[2];
sx q[2];
rz(2.0905154) q[2];
rz(1.4970695) q[3];
sx q[3];
rz(-1.9688508) q[3];
sx q[3];
rz(2.3433949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.368211) q[0];
sx q[0];
rz(-1.0268651) q[0];
sx q[0];
rz(-0.89212242) q[0];
rz(-2.6719773) q[1];
sx q[1];
rz(-0.62925595) q[1];
sx q[1];
rz(-2.3687252) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4887805) q[0];
sx q[0];
rz(-2.1202181) q[0];
sx q[0];
rz(-2.1958952) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8310089) q[2];
sx q[2];
rz(-2.2573244) q[2];
sx q[2];
rz(-1.7169184) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1972496) q[1];
sx q[1];
rz(-2.5641003) q[1];
sx q[1];
rz(-2.5275699) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3350713) q[3];
sx q[3];
rz(-1.3355713) q[3];
sx q[3];
rz(2.1947088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5272687) q[2];
sx q[2];
rz(-1.2246776) q[2];
sx q[2];
rz(0.69250715) q[2];
rz(-2.6912189) q[3];
sx q[3];
rz(-1.9362484) q[3];
sx q[3];
rz(1.6374121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5699128) q[0];
sx q[0];
rz(-2.3501861) q[0];
sx q[0];
rz(0.94863844) q[0];
rz(-1.0665077) q[1];
sx q[1];
rz(-0.98943168) q[1];
sx q[1];
rz(-1.8705286) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9635634) q[0];
sx q[0];
rz(-2.9157713) q[0];
sx q[0];
rz(-1.144335) q[0];
rz(0.022106013) q[2];
sx q[2];
rz(-1.6048631) q[2];
sx q[2];
rz(0.23512041) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3349325) q[1];
sx q[1];
rz(-2.1927823) q[1];
sx q[1];
rz(0.16652624) q[1];
x q[2];
rz(-2.8715114) q[3];
sx q[3];
rz(-2.684786) q[3];
sx q[3];
rz(-2.1169259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5991685) q[2];
sx q[2];
rz(-1.3072661) q[2];
sx q[2];
rz(-0.14031169) q[2];
rz(3.1091651) q[3];
sx q[3];
rz(-2.2166538) q[3];
sx q[3];
rz(1.2062581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5133544) q[0];
sx q[0];
rz(-1.6980549) q[0];
sx q[0];
rz(-2.3640609) q[0];
rz(0.93742049) q[1];
sx q[1];
rz(-1.2317069) q[1];
sx q[1];
rz(0.44313988) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69618689) q[0];
sx q[0];
rz(-1.6517795) q[0];
sx q[0];
rz(0.085025351) q[0];
rz(-pi) q[1];
rz(-0.97340092) q[2];
sx q[2];
rz(-1.8466766) q[2];
sx q[2];
rz(1.2593002) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.41105553) q[1];
sx q[1];
rz(-1.2403508) q[1];
sx q[1];
rz(2.4871268) q[1];
rz(-2.0432161) q[3];
sx q[3];
rz(-2.7948423) q[3];
sx q[3];
rz(-0.43671331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0605269) q[2];
sx q[2];
rz(-2.492283) q[2];
sx q[2];
rz(-0.12760663) q[2];
rz(0.29868948) q[3];
sx q[3];
rz(-1.1726215) q[3];
sx q[3];
rz(2.7711788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5304607) q[0];
sx q[0];
rz(-1.8831384) q[0];
sx q[0];
rz(0.25181121) q[0];
rz(-0.61027377) q[1];
sx q[1];
rz(-0.88519575) q[1];
sx q[1];
rz(-2.0559678) q[1];
rz(-3.0771607) q[2];
sx q[2];
rz(-2.2885135) q[2];
sx q[2];
rz(2.0533333) q[2];
rz(-1.1202624) q[3];
sx q[3];
rz(-1.7623175) q[3];
sx q[3];
rz(0.13468066) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
