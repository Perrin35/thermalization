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
rz(0.62975878) q[0];
sx q[0];
rz(3.44343) q[0];
sx q[0];
rz(11.259196) q[0];
rz(-3.6699927) q[1];
sx q[1];
rz(3.9730605) q[1];
sx q[1];
rz(12.107036) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2506901) q[0];
sx q[0];
rz(-2.8807682) q[0];
sx q[0];
rz(-1.7623259) q[0];
x q[1];
rz(-1.8474962) q[2];
sx q[2];
rz(-1.7613698) q[2];
sx q[2];
rz(-1.1016359) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3078008) q[1];
sx q[1];
rz(-1.7437309) q[1];
sx q[1];
rz(0.61638919) q[1];
x q[2];
rz(1.5993871) q[3];
sx q[3];
rz(-0.66900476) q[3];
sx q[3];
rz(-2.1137848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1186195) q[2];
sx q[2];
rz(-1.8034673) q[2];
sx q[2];
rz(-2.5845134) q[2];
rz(-0.45595512) q[3];
sx q[3];
rz(-2.3796701) q[3];
sx q[3];
rz(-0.60321155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2073681) q[0];
sx q[0];
rz(-0.71894431) q[0];
sx q[0];
rz(1.7508605) q[0];
rz(-1.8234183) q[1];
sx q[1];
rz(-1.428182) q[1];
sx q[1];
rz(-0.21600977) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7129546) q[0];
sx q[0];
rz(-0.99073505) q[0];
sx q[0];
rz(2.7614276) q[0];
rz(0.26311263) q[2];
sx q[2];
rz(-2.7091658) q[2];
sx q[2];
rz(-0.95460923) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60571721) q[1];
sx q[1];
rz(-0.89369666) q[1];
sx q[1];
rz(0.89800055) q[1];
rz(-1.9253821) q[3];
sx q[3];
rz(-2.5719593) q[3];
sx q[3];
rz(-0.58651741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7184489) q[2];
sx q[2];
rz(-2.7687912) q[2];
sx q[2];
rz(0.049840363) q[2];
rz(2.5981264) q[3];
sx q[3];
rz(-1.1246357) q[3];
sx q[3];
rz(-2.0487093) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8922888) q[0];
sx q[0];
rz(-2.0396621) q[0];
sx q[0];
rz(-2.1732543) q[0];
rz(3.1210506) q[1];
sx q[1];
rz(-1.7612709) q[1];
sx q[1];
rz(-1.7113908) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1002625) q[0];
sx q[0];
rz(-1.593129) q[0];
sx q[0];
rz(0.64133994) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5495785) q[2];
sx q[2];
rz(-0.63767728) q[2];
sx q[2];
rz(-2.9108832) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.80950981) q[1];
sx q[1];
rz(-0.44594279) q[1];
sx q[1];
rz(2.6757702) q[1];
rz(-pi) q[2];
rz(2.7647721) q[3];
sx q[3];
rz(-1.539789) q[3];
sx q[3];
rz(3.0135807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.7057544) q[2];
sx q[2];
rz(-0.64283723) q[2];
sx q[2];
rz(1.7819116) q[2];
rz(-2.6333574) q[3];
sx q[3];
rz(-1.619092) q[3];
sx q[3];
rz(1.9967509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(1.085676) q[0];
sx q[0];
rz(-0.29656947) q[0];
sx q[0];
rz(-1.0067518) q[0];
rz(0.83167568) q[1];
sx q[1];
rz(-0.74200231) q[1];
sx q[1];
rz(1.1078018) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2058973) q[0];
sx q[0];
rz(-0.95608053) q[0];
sx q[0];
rz(2.6877235) q[0];
rz(2.7659589) q[2];
sx q[2];
rz(-1.8058597) q[2];
sx q[2];
rz(1.7441074) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.133114) q[1];
sx q[1];
rz(-1.8777056) q[1];
sx q[1];
rz(0.020773356) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6964968) q[3];
sx q[3];
rz(-1.8254878) q[3];
sx q[3];
rz(-1.2459918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45336777) q[2];
sx q[2];
rz(-2.0757165) q[2];
sx q[2];
rz(0.31309703) q[2];
rz(0.39737663) q[3];
sx q[3];
rz(-1.671096) q[3];
sx q[3];
rz(-2.9562922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6322286) q[0];
sx q[0];
rz(-2.2780184) q[0];
sx q[0];
rz(-2.2160227) q[0];
rz(-2.6717692) q[1];
sx q[1];
rz(-1.3146105) q[1];
sx q[1];
rz(2.125461) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6189656) q[0];
sx q[0];
rz(-2.5182808) q[0];
sx q[0];
rz(-0.81617426) q[0];
x q[1];
rz(-2.2158898) q[2];
sx q[2];
rz(-2.7113554) q[2];
sx q[2];
rz(-1.0412585) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.092246902) q[1];
sx q[1];
rz(-1.3692442) q[1];
sx q[1];
rz(-1.818403) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8941325) q[3];
sx q[3];
rz(-1.5405021) q[3];
sx q[3];
rz(0.39112511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65285811) q[2];
sx q[2];
rz(-0.32565871) q[2];
sx q[2];
rz(2.6540836) q[2];
rz(-2.2440535) q[3];
sx q[3];
rz(-1.2582015) q[3];
sx q[3];
rz(-2.0770309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2010736) q[0];
sx q[0];
rz(-0.93795332) q[0];
sx q[0];
rz(2.4679389) q[0];
rz(1.9081839) q[1];
sx q[1];
rz(-2.1234832) q[1];
sx q[1];
rz(2.3755551) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4318643) q[0];
sx q[0];
rz(-2.0805295) q[0];
sx q[0];
rz(-2.8317189) q[0];
rz(0.17579097) q[2];
sx q[2];
rz(-1.5002709) q[2];
sx q[2];
rz(2.8688237) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.627558) q[1];
sx q[1];
rz(-2.6366391) q[1];
sx q[1];
rz(1.2673753) q[1];
rz(-pi) q[2];
rz(0.30575846) q[3];
sx q[3];
rz(-1.7019203) q[3];
sx q[3];
rz(0.46520376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.172714) q[2];
sx q[2];
rz(-1.4977027) q[2];
sx q[2];
rz(-2.3806351) q[2];
rz(2.739665) q[3];
sx q[3];
rz(-2.8765078) q[3];
sx q[3];
rz(-2.5529805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6589979) q[0];
sx q[0];
rz(-0.75356475) q[0];
sx q[0];
rz(1.448451) q[0];
rz(0.50085577) q[1];
sx q[1];
rz(-1.294699) q[1];
sx q[1];
rz(-2.3282611) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7194448) q[0];
sx q[0];
rz(-1.1851743) q[0];
sx q[0];
rz(-1.411012) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3766889) q[2];
sx q[2];
rz(-0.55667215) q[2];
sx q[2];
rz(-1.0355606) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.2351027) q[1];
sx q[1];
rz(-1.1336599) q[1];
sx q[1];
rz(-1.7227931) q[1];
x q[2];
rz(1.1239979) q[3];
sx q[3];
rz(-0.13544336) q[3];
sx q[3];
rz(-2.9970084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.27998754) q[2];
sx q[2];
rz(-0.6051175) q[2];
sx q[2];
rz(0.32312265) q[2];
rz(1.4019639) q[3];
sx q[3];
rz(-1.2819382) q[3];
sx q[3];
rz(1.7590947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0793656) q[0];
sx q[0];
rz(-1.0492188) q[0];
sx q[0];
rz(-2.3754689) q[0];
rz(0.8194204) q[1];
sx q[1];
rz(-1.0304281) q[1];
sx q[1];
rz(-1.5310418) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3643506) q[0];
sx q[0];
rz(-1.3421699) q[0];
sx q[0];
rz(1.6407439) q[0];
x q[1];
rz(0.28240164) q[2];
sx q[2];
rz(-2.4666383) q[2];
sx q[2];
rz(-2.9815428) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4041679) q[1];
sx q[1];
rz(-2.1497207) q[1];
sx q[1];
rz(-2.3217724) q[1];
x q[2];
rz(1.895745) q[3];
sx q[3];
rz(-1.2940931) q[3];
sx q[3];
rz(0.70223728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7649585) q[2];
sx q[2];
rz(-2.9464293) q[2];
sx q[2];
rz(-2.7435319) q[2];
rz(-2.5352488) q[3];
sx q[3];
rz(-1.5747993) q[3];
sx q[3];
rz(-0.91355598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9048318) q[0];
sx q[0];
rz(-0.30192152) q[0];
sx q[0];
rz(2.6171369) q[0];
rz(2.4728788) q[1];
sx q[1];
rz(-1.6122183) q[1];
sx q[1];
rz(-1.761577) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3052914) q[0];
sx q[0];
rz(-2.143465) q[0];
sx q[0];
rz(0.051866626) q[0];
rz(-pi) q[1];
rz(-0.30461208) q[2];
sx q[2];
rz(-1.1122983) q[2];
sx q[2];
rz(0.079041399) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6898247) q[1];
sx q[1];
rz(-2.5344026) q[1];
sx q[1];
rz(-2.3535806) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1942176) q[3];
sx q[3];
rz(-2.0717952) q[3];
sx q[3];
rz(0.11160103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9809208) q[2];
sx q[2];
rz(-2.394684) q[2];
sx q[2];
rz(2.6778632) q[2];
rz(-1.7175698) q[3];
sx q[3];
rz(-1.0878891) q[3];
sx q[3];
rz(-0.033528479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56228191) q[0];
sx q[0];
rz(-0.71826851) q[0];
sx q[0];
rz(-0.41807362) q[0];
rz(2.383291) q[1];
sx q[1];
rz(-2.1577991) q[1];
sx q[1];
rz(-1.2850579) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3089963) q[0];
sx q[0];
rz(-0.65387883) q[0];
sx q[0];
rz(0.12419463) q[0];
rz(-pi) q[1];
rz(1.1056771) q[2];
sx q[2];
rz(-2.4659202) q[2];
sx q[2];
rz(0.0059520324) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.27800769) q[1];
sx q[1];
rz(-1.3787495) q[1];
sx q[1];
rz(2.5878367) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0009406) q[3];
sx q[3];
rz(-0.39467282) q[3];
sx q[3];
rz(-0.81854023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9717504) q[2];
sx q[2];
rz(-1.0047793) q[2];
sx q[2];
rz(-0.24610914) q[2];
rz(-1.2717815) q[3];
sx q[3];
rz(-1.566794) q[3];
sx q[3];
rz(1.3625712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1590189) q[0];
sx q[0];
rz(-1.6652501) q[0];
sx q[0];
rz(-0.2022947) q[0];
rz(-0.62494878) q[1];
sx q[1];
rz(-2.2849871) q[1];
sx q[1];
rz(-2.7402592) q[1];
rz(2.305793) q[2];
sx q[2];
rz(-0.27793548) q[2];
sx q[2];
rz(-0.64151573) q[2];
rz(-2.940964) q[3];
sx q[3];
rz(-1.2617785) q[3];
sx q[3];
rz(-2.3599335) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
