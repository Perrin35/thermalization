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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.280312) q[0];
sx q[0];
rz(-0.93398184) q[0];
sx q[0];
rz(-2.0122819) q[0];
rz(-pi) q[1];
rz(0.91602625) q[2];
sx q[2];
rz(-1.6189304) q[2];
sx q[2];
rz(-2.2544238) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4976884) q[1];
sx q[1];
rz(-2.4329281) q[1];
sx q[1];
rz(-0.24073118) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2774062) q[3];
sx q[3];
rz(-1.7134616) q[3];
sx q[3];
rz(-2.8418926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76838905) q[2];
sx q[2];
rz(-1.6037805) q[2];
sx q[2];
rz(-2.0453889) q[2];
rz(1.0783892) q[3];
sx q[3];
rz(-0.53467852) q[3];
sx q[3];
rz(2.6064742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.815149) q[0];
sx q[0];
rz(-2.09477) q[0];
sx q[0];
rz(-0.4775508) q[0];
rz(2.1304456) q[1];
sx q[1];
rz(-0.54713455) q[1];
sx q[1];
rz(-2.0571041) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3318429) q[0];
sx q[0];
rz(-1.6879775) q[0];
sx q[0];
rz(-1.2007037) q[0];
rz(-pi) q[1];
rz(-1.162503) q[2];
sx q[2];
rz(-1.5973685) q[2];
sx q[2];
rz(2.0808737) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0140927) q[1];
sx q[1];
rz(-1.0261826) q[1];
sx q[1];
rz(1.3633534) q[1];
rz(-pi) q[2];
rz(-1.7803187) q[3];
sx q[3];
rz(-2.2993436) q[3];
sx q[3];
rz(1.9467348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5227205) q[2];
sx q[2];
rz(-2.7174157) q[2];
sx q[2];
rz(-2.0514533) q[2];
rz(2.5777396) q[3];
sx q[3];
rz(-0.68958759) q[3];
sx q[3];
rz(-3.1088945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3801333) q[0];
sx q[0];
rz(-3.1205803) q[0];
sx q[0];
rz(1.6047961) q[0];
rz(1.3179294) q[1];
sx q[1];
rz(-2.047796) q[1];
sx q[1];
rz(0.39753786) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4039392) q[0];
sx q[0];
rz(-1.353147) q[0];
sx q[0];
rz(0.91380878) q[0];
rz(2.1486542) q[2];
sx q[2];
rz(-1.8210337) q[2];
sx q[2];
rz(-1.9577574) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.51008114) q[1];
sx q[1];
rz(-2.1891174) q[1];
sx q[1];
rz(-0.31429283) q[1];
rz(2.0769172) q[3];
sx q[3];
rz(-1.1017513) q[3];
sx q[3];
rz(2.7739547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90618187) q[2];
sx q[2];
rz(-1.9637039) q[2];
sx q[2];
rz(0.57514352) q[2];
rz(2.4681674) q[3];
sx q[3];
rz(-1.5410825) q[3];
sx q[3];
rz(1.5448236) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49066082) q[0];
sx q[0];
rz(-1.1216102) q[0];
sx q[0];
rz(0.90890539) q[0];
rz(3.124369) q[1];
sx q[1];
rz(-2.6135018) q[1];
sx q[1];
rz(-0.012103279) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72649137) q[0];
sx q[0];
rz(-2.4612777) q[0];
sx q[0];
rz(2.6136287) q[0];
x q[1];
rz(-2.1075756) q[2];
sx q[2];
rz(-0.99269789) q[2];
sx q[2];
rz(-0.54852099) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3084532) q[1];
sx q[1];
rz(-2.1670682) q[1];
sx q[1];
rz(-2.8677529) q[1];
rz(-3.0148325) q[3];
sx q[3];
rz(-1.5258299) q[3];
sx q[3];
rz(0.25925207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0688613) q[2];
sx q[2];
rz(-2.8673745) q[2];
sx q[2];
rz(-2.7552628) q[2];
rz(-0.38935152) q[3];
sx q[3];
rz(-1.5334305) q[3];
sx q[3];
rz(2.1336335) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7858081) q[0];
sx q[0];
rz(-2.4171827) q[0];
sx q[0];
rz(0.32387787) q[0];
rz(-2.7549699) q[1];
sx q[1];
rz(-1.7522248) q[1];
sx q[1];
rz(-1.5868384) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1607035) q[0];
sx q[0];
rz(-1.6352773) q[0];
sx q[0];
rz(2.4095201) q[0];
x q[1];
rz(-2.1542495) q[2];
sx q[2];
rz(-0.61827055) q[2];
sx q[2];
rz(1.5702919) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7534947) q[1];
sx q[1];
rz(-0.83067229) q[1];
sx q[1];
rz(0.64708935) q[1];
rz(2.0915786) q[3];
sx q[3];
rz(-0.84283295) q[3];
sx q[3];
rz(-1.4340374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3891478) q[2];
sx q[2];
rz(-0.84263313) q[2];
sx q[2];
rz(-3.0774934) q[2];
rz(-0.60424232) q[3];
sx q[3];
rz(-1.3925545) q[3];
sx q[3];
rz(-1.909168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9391249) q[0];
sx q[0];
rz(-0.61102837) q[0];
sx q[0];
rz(-0.87919277) q[0];
rz(-2.7912256) q[1];
sx q[1];
rz(-1.8354974) q[1];
sx q[1];
rz(2.9894357) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2911173) q[0];
sx q[0];
rz(-0.88248173) q[0];
sx q[0];
rz(1.2173714) q[0];
rz(1.2142356) q[2];
sx q[2];
rz(-2.6313553) q[2];
sx q[2];
rz(-1.7348926) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.31649905) q[1];
sx q[1];
rz(-0.48227019) q[1];
sx q[1];
rz(2.6939166) q[1];
rz(2.426126) q[3];
sx q[3];
rz(-2.5757534) q[3];
sx q[3];
rz(-2.6239606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0396314) q[2];
sx q[2];
rz(-0.32932082) q[2];
sx q[2];
rz(-2.9015818) q[2];
rz(-2.4890684) q[3];
sx q[3];
rz(-1.1381166) q[3];
sx q[3];
rz(-1.8664546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0535102) q[0];
sx q[0];
rz(-1.2367915) q[0];
sx q[0];
rz(0.85813338) q[0];
rz(1.3872604) q[1];
sx q[1];
rz(-0.45448449) q[1];
sx q[1];
rz(2.8087356) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3368206) q[0];
sx q[0];
rz(-2.4737278) q[0];
sx q[0];
rz(-1.0475545) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9152739) q[2];
sx q[2];
rz(-0.25626015) q[2];
sx q[2];
rz(1.8654491) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7499125) q[1];
sx q[1];
rz(-1.1521764) q[1];
sx q[1];
rz(1.7961572) q[1];
x q[2];
rz(2.4868785) q[3];
sx q[3];
rz(-1.2867478) q[3];
sx q[3];
rz(-2.9971307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.7994999) q[2];
sx q[2];
rz(-2.9436593) q[2];
sx q[2];
rz(-2.0905154) q[2];
rz(-1.6445232) q[3];
sx q[3];
rz(-1.1727419) q[3];
sx q[3];
rz(0.79819775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.368211) q[0];
sx q[0];
rz(-2.1147275) q[0];
sx q[0];
rz(-0.89212242) q[0];
rz(-2.6719773) q[1];
sx q[1];
rz(-2.5123367) q[1];
sx q[1];
rz(2.3687252) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55760705) q[0];
sx q[0];
rz(-2.0933042) q[0];
sx q[0];
rz(2.4948392) q[0];
rz(-1.8310089) q[2];
sx q[2];
rz(-0.88426829) q[2];
sx q[2];
rz(1.4246743) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6438021) q[1];
sx q[1];
rz(-2.0333159) q[1];
sx q[1];
rz(-1.9299354) q[1];
rz(-pi) q[2];
rz(1.8065213) q[3];
sx q[3];
rz(-1.8060214) q[3];
sx q[3];
rz(-0.94688382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5272687) q[2];
sx q[2];
rz(-1.9169151) q[2];
sx q[2];
rz(0.69250715) q[2];
rz(0.45037371) q[3];
sx q[3];
rz(-1.2053442) q[3];
sx q[3];
rz(1.5041806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57167989) q[0];
sx q[0];
rz(-2.3501861) q[0];
sx q[0];
rz(0.94863844) q[0];
rz(-1.0665077) q[1];
sx q[1];
rz(-2.152161) q[1];
sx q[1];
rz(1.8705286) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5273351) q[0];
sx q[0];
rz(-1.7760881) q[0];
sx q[0];
rz(0.094747748) q[0];
x q[1];
rz(0.99540751) q[2];
sx q[2];
rz(-0.040608309) q[2];
sx q[2];
rz(2.3307068) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8066602) q[1];
sx q[1];
rz(-2.1927823) q[1];
sx q[1];
rz(-0.16652624) q[1];
rz(0.27008121) q[3];
sx q[3];
rz(-2.684786) q[3];
sx q[3];
rz(-2.1169259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5991685) q[2];
sx q[2];
rz(-1.3072661) q[2];
sx q[2];
rz(-0.14031169) q[2];
rz(3.1091651) q[3];
sx q[3];
rz(-2.2166538) q[3];
sx q[3];
rz(-1.9353346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62823826) q[0];
sx q[0];
rz(-1.4435377) q[0];
sx q[0];
rz(0.7775318) q[0];
rz(2.2041722) q[1];
sx q[1];
rz(-1.2317069) q[1];
sx q[1];
rz(-0.44313988) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4454058) q[0];
sx q[0];
rz(-1.6517795) q[0];
sx q[0];
rz(3.0565673) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1681917) q[2];
sx q[2];
rz(-1.8466766) q[2];
sx q[2];
rz(1.8822924) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7305371) q[1];
sx q[1];
rz(-1.9012419) q[1];
sx q[1];
rz(2.4871268) q[1];
x q[2];
rz(2.9786213) q[3];
sx q[3];
rz(-1.8782369) q[3];
sx q[3];
rz(-2.2071537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0810658) q[2];
sx q[2];
rz(-0.64930969) q[2];
sx q[2];
rz(-0.12760663) q[2];
rz(-0.29868948) q[3];
sx q[3];
rz(-1.9689711) q[3];
sx q[3];
rz(2.7711788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5304607) q[0];
sx q[0];
rz(-1.2584542) q[0];
sx q[0];
rz(-2.8897814) q[0];
rz(0.61027377) q[1];
sx q[1];
rz(-2.2563969) q[1];
sx q[1];
rz(1.0856249) q[1];
rz(2.2895428) q[2];
sx q[2];
rz(-1.5222737) q[2];
sx q[2];
rz(0.52494502) q[2];
rz(-2.0213303) q[3];
sx q[3];
rz(-1.3792752) q[3];
sx q[3];
rz(-3.006912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
