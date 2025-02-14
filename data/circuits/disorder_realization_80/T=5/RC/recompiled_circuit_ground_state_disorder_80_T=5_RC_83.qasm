OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4898981) q[0];
sx q[0];
rz(-2.2597921) q[0];
sx q[0];
rz(-2.9969126) q[0];
rz(-2.7235003) q[1];
sx q[1];
rz(-3.0378208) q[1];
sx q[1];
rz(1.6296847) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1967752) q[0];
sx q[0];
rz(-0.51048764) q[0];
sx q[0];
rz(1.8148242) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6250526) q[2];
sx q[2];
rz(-2.3035604) q[2];
sx q[2];
rz(2.4861479) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.92951951) q[1];
sx q[1];
rz(-1.9572502) q[1];
sx q[1];
rz(0.35915908) q[1];
rz(-pi) q[2];
rz(1.1505262) q[3];
sx q[3];
rz(-1.0967249) q[3];
sx q[3];
rz(-1.898511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.31643733) q[2];
sx q[2];
rz(-1.8708159) q[2];
sx q[2];
rz(-2.4862508) q[2];
rz(1.3842899) q[3];
sx q[3];
rz(-2.1770848) q[3];
sx q[3];
rz(-2.3206319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4493988) q[0];
sx q[0];
rz(-1.3757502) q[0];
sx q[0];
rz(-0.50814116) q[0];
rz(1.3805768) q[1];
sx q[1];
rz(-0.5813798) q[1];
sx q[1];
rz(1.2095721) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.094645) q[0];
sx q[0];
rz(-0.98804615) q[0];
sx q[0];
rz(-0.21999448) q[0];
x q[1];
rz(-0.31866535) q[2];
sx q[2];
rz(-0.056431596) q[2];
sx q[2];
rz(2.3669764) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.74871127) q[1];
sx q[1];
rz(-1.4429174) q[1];
sx q[1];
rz(0.46462469) q[1];
rz(-pi) q[2];
rz(1.0367514) q[3];
sx q[3];
rz(-2.4656714) q[3];
sx q[3];
rz(0.40646857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2035344) q[2];
sx q[2];
rz(-1.7455696) q[2];
sx q[2];
rz(-0.61689845) q[2];
rz(1.5361702) q[3];
sx q[3];
rz(-2.0974396) q[3];
sx q[3];
rz(2.6507586) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9817552) q[0];
sx q[0];
rz(-1.497739) q[0];
sx q[0];
rz(0.72506654) q[0];
rz(-1.8720576) q[1];
sx q[1];
rz(-0.67765647) q[1];
sx q[1];
rz(-2.1824172) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3676783) q[0];
sx q[0];
rz(-1.1836021) q[0];
sx q[0];
rz(-2.0810803) q[0];
rz(-pi) q[1];
rz(1.1821907) q[2];
sx q[2];
rz(-1.1827785) q[2];
sx q[2];
rz(-1.4997375) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37466418) q[1];
sx q[1];
rz(-2.9247724) q[1];
sx q[1];
rz(-1.0105074) q[1];
rz(-pi) q[2];
x q[2];
rz(2.972288) q[3];
sx q[3];
rz(-1.7506071) q[3];
sx q[3];
rz(1.394751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3265257) q[2];
sx q[2];
rz(-1.2667789) q[2];
sx q[2];
rz(2.8845924) q[2];
rz(-2.0604996) q[3];
sx q[3];
rz(-2.6205781) q[3];
sx q[3];
rz(1.5787554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0839888) q[0];
sx q[0];
rz(-0.27844089) q[0];
sx q[0];
rz(1.1778911) q[0];
rz(0.15011694) q[1];
sx q[1];
rz(-2.6094486) q[1];
sx q[1];
rz(-1.5163126) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35203194) q[0];
sx q[0];
rz(-1.9012678) q[0];
sx q[0];
rz(0.29833692) q[0];
rz(-pi) q[1];
rz(-0.37090918) q[2];
sx q[2];
rz(-0.87693767) q[2];
sx q[2];
rz(-2.3313076) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8575688) q[1];
sx q[1];
rz(-0.34115667) q[1];
sx q[1];
rz(2.8400303) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7606973) q[3];
sx q[3];
rz(-0.93578029) q[3];
sx q[3];
rz(0.39228016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2991422) q[2];
sx q[2];
rz(-0.56988847) q[2];
sx q[2];
rz(-2.2184856) q[2];
rz(-2.182377) q[3];
sx q[3];
rz(-1.0126637) q[3];
sx q[3];
rz(-0.21624163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32090309) q[0];
sx q[0];
rz(-0.90383363) q[0];
sx q[0];
rz(2.2586816) q[0];
rz(-1.8461022) q[1];
sx q[1];
rz(-2.3675282) q[1];
sx q[1];
rz(0.49497089) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6225016) q[0];
sx q[0];
rz(-1.8880196) q[0];
sx q[0];
rz(-0.30479161) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70303838) q[2];
sx q[2];
rz(-1.993317) q[2];
sx q[2];
rz(-0.51681821) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.40303142) q[1];
sx q[1];
rz(-0.90442362) q[1];
sx q[1];
rz(-1.6825466) q[1];
x q[2];
rz(0.35843973) q[3];
sx q[3];
rz(-1.0762941) q[3];
sx q[3];
rz(0.93264893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3317269) q[2];
sx q[2];
rz(-0.11955424) q[2];
sx q[2];
rz(1.8734107) q[2];
rz(-0.082503334) q[3];
sx q[3];
rz(-1.5039597) q[3];
sx q[3];
rz(-0.65703195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4853915) q[0];
sx q[0];
rz(-2.9347561) q[0];
sx q[0];
rz(1.50151) q[0];
rz(1.062475) q[1];
sx q[1];
rz(-1.1524009) q[1];
sx q[1];
rz(-1.9662439) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51477434) q[0];
sx q[0];
rz(-1.2668224) q[0];
sx q[0];
rz(-1.8016812) q[0];
rz(2.4480216) q[2];
sx q[2];
rz(-1.3422924) q[2];
sx q[2];
rz(-2.9442043) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3835241) q[1];
sx q[1];
rz(-1.822377) q[1];
sx q[1];
rz(-1.4497767) q[1];
rz(-pi) q[2];
rz(1.2731984) q[3];
sx q[3];
rz(-0.83286017) q[3];
sx q[3];
rz(-2.0119865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33092734) q[2];
sx q[2];
rz(-1.8570447) q[2];
sx q[2];
rz(0.98768273) q[2];
rz(0.87743131) q[3];
sx q[3];
rz(-2.8743447) q[3];
sx q[3];
rz(0.73006829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.1838609) q[0];
sx q[0];
rz(-1.5393625) q[0];
sx q[0];
rz(-3.0027332) q[0];
rz(1.9271556) q[1];
sx q[1];
rz(-1.5077533) q[1];
sx q[1];
rz(-2.3249998) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2557295) q[0];
sx q[0];
rz(-0.22561377) q[0];
sx q[0];
rz(-1.4855235) q[0];
rz(-1.310552) q[2];
sx q[2];
rz(-2.3219371) q[2];
sx q[2];
rz(-1.9323088) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5467756) q[1];
sx q[1];
rz(-2.4218514) q[1];
sx q[1];
rz(-1.627658) q[1];
x q[2];
rz(0.84202977) q[3];
sx q[3];
rz(-2.2737164) q[3];
sx q[3];
rz(1.3852392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.82178086) q[2];
sx q[2];
rz(-0.37415162) q[2];
sx q[2];
rz(-2.6569341) q[2];
rz(-2.7759077) q[3];
sx q[3];
rz(-1.7771143) q[3];
sx q[3];
rz(1.1588089) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1300238) q[0];
sx q[0];
rz(-0.31112177) q[0];
sx q[0];
rz(3.044361) q[0];
rz(-0.10487996) q[1];
sx q[1];
rz(-2.361894) q[1];
sx q[1];
rz(-1.3146776) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.437156) q[0];
sx q[0];
rz(-1.5644685) q[0];
sx q[0];
rz(-3.1392211) q[0];
rz(-pi) q[1];
rz(0.22366053) q[2];
sx q[2];
rz(-1.3350492) q[2];
sx q[2];
rz(-3.0652114) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2995616) q[1];
sx q[1];
rz(-0.52134575) q[1];
sx q[1];
rz(-0.10592769) q[1];
x q[2];
rz(1.0339853) q[3];
sx q[3];
rz(-1.4958463) q[3];
sx q[3];
rz(-1.0190462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17681992) q[2];
sx q[2];
rz(-1.3975881) q[2];
sx q[2];
rz(-0.1532661) q[2];
rz(-1.2166474) q[3];
sx q[3];
rz(-2.9235268) q[3];
sx q[3];
rz(0.6161859) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1185054) q[0];
sx q[0];
rz(-3.1361134) q[0];
sx q[0];
rz(1.5047005) q[0];
rz(2.3432689) q[1];
sx q[1];
rz(-1.6183805) q[1];
sx q[1];
rz(-2.8062779) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.766205) q[0];
sx q[0];
rz(-1.6339301) q[0];
sx q[0];
rz(-1.8324018) q[0];
rz(1.5463355) q[2];
sx q[2];
rz(-0.82055295) q[2];
sx q[2];
rz(-2.1945107) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1095417) q[1];
sx q[1];
rz(-0.45278835) q[1];
sx q[1];
rz(0.70346634) q[1];
rz(-pi) q[2];
rz(-2.1481403) q[3];
sx q[3];
rz(-1.7029188) q[3];
sx q[3];
rz(-2.7293049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.12751427) q[2];
sx q[2];
rz(-1.9731584) q[2];
sx q[2];
rz(0.43759313) q[2];
rz(2.0994999) q[3];
sx q[3];
rz(-1.7362678) q[3];
sx q[3];
rz(-0.54245814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9893148) q[0];
sx q[0];
rz(-0.88459009) q[0];
sx q[0];
rz(-0.44961318) q[0];
rz(2.1647029) q[1];
sx q[1];
rz(-1.6376817) q[1];
sx q[1];
rz(-0.50122112) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4924663) q[0];
sx q[0];
rz(-1.7607795) q[0];
sx q[0];
rz(-1.282592) q[0];
x q[1];
rz(2.0222763) q[2];
sx q[2];
rz(-2.6767133) q[2];
sx q[2];
rz(1.3550188) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.93443643) q[1];
sx q[1];
rz(-2.0324576) q[1];
sx q[1];
rz(0.2646048) q[1];
rz(-pi) q[2];
rz(-2.7990255) q[3];
sx q[3];
rz(-1.1293355) q[3];
sx q[3];
rz(1.6395237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1180798) q[2];
sx q[2];
rz(-2.2492275) q[2];
sx q[2];
rz(0.21772131) q[2];
rz(1.2348385) q[3];
sx q[3];
rz(-0.66399884) q[3];
sx q[3];
rz(-2.2209404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6476743) q[0];
sx q[0];
rz(-1.4287345) q[0];
sx q[0];
rz(-2.7609974) q[0];
rz(-2.4299798) q[1];
sx q[1];
rz(-1.8743534) q[1];
sx q[1];
rz(1.4030917) q[1];
rz(-2.6720302) q[2];
sx q[2];
rz(-1.6161412) q[2];
sx q[2];
rz(-1.8117088) q[2];
rz(-2.7357581) q[3];
sx q[3];
rz(-1.3530227) q[3];
sx q[3];
rz(1.4128662) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
