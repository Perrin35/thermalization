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
rz(-2.45911) q[0];
sx q[0];
rz(-0.64829666) q[0];
sx q[0];
rz(0.034870235) q[0];
rz(-1.7416481) q[1];
sx q[1];
rz(-2.8748684) q[1];
sx q[1];
rz(1.9884225) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1138251) q[0];
sx q[0];
rz(-1.5875289) q[0];
sx q[0];
rz(-1.2196779) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2049293) q[2];
sx q[2];
rz(-1.7589169) q[2];
sx q[2];
rz(0.29292303) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0583078) q[1];
sx q[1];
rz(-2.1911006) q[1];
sx q[1];
rz(-2.9172667) q[1];
rz(3.127159) q[3];
sx q[3];
rz(-2.0216515) q[3];
sx q[3];
rz(0.23896398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1109041) q[2];
sx q[2];
rz(-1.2653753) q[2];
sx q[2];
rz(0.14070025) q[2];
rz(-1.5594907) q[3];
sx q[3];
rz(-2.9086106) q[3];
sx q[3];
rz(-1.8039907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88386184) q[0];
sx q[0];
rz(-0.95153874) q[0];
sx q[0];
rz(2.0146433) q[0];
rz(-1.1832712) q[1];
sx q[1];
rz(-1.1401581) q[1];
sx q[1];
rz(0.3068876) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81079262) q[0];
sx q[0];
rz(-1.188827) q[0];
sx q[0];
rz(1.6195787) q[0];
rz(-pi) q[1];
rz(0.82123001) q[2];
sx q[2];
rz(-2.3013287) q[2];
sx q[2];
rz(2.822838) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3813016) q[1];
sx q[1];
rz(-2.8861001) q[1];
sx q[1];
rz(-2.5653272) q[1];
x q[2];
rz(1.0829665) q[3];
sx q[3];
rz(-1.6681507) q[3];
sx q[3];
rz(3.0157582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41944501) q[2];
sx q[2];
rz(-0.38997969) q[2];
sx q[2];
rz(-1.5004213) q[2];
rz(1.8309343) q[3];
sx q[3];
rz(-2.1597517) q[3];
sx q[3];
rz(2.380044) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0797043) q[0];
sx q[0];
rz(-2.5766928) q[0];
sx q[0];
rz(1.0796219) q[0];
rz(2.9234746) q[1];
sx q[1];
rz(-1.7354542) q[1];
sx q[1];
rz(-1.1880147) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2750888) q[0];
sx q[0];
rz(-1.554477) q[0];
sx q[0];
rz(-0.62243502) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57680371) q[2];
sx q[2];
rz(-1.7110023) q[2];
sx q[2];
rz(1.1568251) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0674141) q[1];
sx q[1];
rz(-1.4715096) q[1];
sx q[1];
rz(-1.5669785) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5651592) q[3];
sx q[3];
rz(-1.7923017) q[3];
sx q[3];
rz(3.053772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5334566) q[2];
sx q[2];
rz(-0.44508219) q[2];
sx q[2];
rz(0.76370561) q[2];
rz(0.25460994) q[3];
sx q[3];
rz(-2.2201241) q[3];
sx q[3];
rz(0.34847954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3449645) q[0];
sx q[0];
rz(-0.94709713) q[0];
sx q[0];
rz(-3.1297041) q[0];
rz(-2.1628926) q[1];
sx q[1];
rz(-1.782676) q[1];
sx q[1];
rz(2.0174644) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55292144) q[0];
sx q[0];
rz(-1.8055989) q[0];
sx q[0];
rz(2.7166883) q[0];
x q[1];
rz(2.5225183) q[2];
sx q[2];
rz(-0.67712578) q[2];
sx q[2];
rz(-0.4437333) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.970579) q[1];
sx q[1];
rz(-0.56736028) q[1];
sx q[1];
rz(2.2544276) q[1];
rz(-2.4120578) q[3];
sx q[3];
rz(-1.9771061) q[3];
sx q[3];
rz(-1.1247579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0594788) q[2];
sx q[2];
rz(-0.79136005) q[2];
sx q[2];
rz(2.3477083) q[2];
rz(-1.3715028) q[3];
sx q[3];
rz(-0.81727782) q[3];
sx q[3];
rz(2.1577458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3802948) q[0];
sx q[0];
rz(-1.607432) q[0];
sx q[0];
rz(-2.5979331) q[0];
rz(-1.1477973) q[1];
sx q[1];
rz(-2.5416083) q[1];
sx q[1];
rz(2.0701087) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54052559) q[0];
sx q[0];
rz(-0.50583848) q[0];
sx q[0];
rz(0.95212014) q[0];
rz(-pi) q[1];
rz(0.29419326) q[2];
sx q[2];
rz(-1.4227267) q[2];
sx q[2];
rz(-1.3180871) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8036441) q[1];
sx q[1];
rz(-2.0211086) q[1];
sx q[1];
rz(1.2195576) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50258639) q[3];
sx q[3];
rz(-2.2690176) q[3];
sx q[3];
rz(1.0909441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.078309623) q[2];
sx q[2];
rz(-2.6368243) q[2];
sx q[2];
rz(0.90739352) q[2];
rz(-3.1388969) q[3];
sx q[3];
rz(-1.7406274) q[3];
sx q[3];
rz(1.2292954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9232848) q[0];
sx q[0];
rz(-1.8508428) q[0];
sx q[0];
rz(1.1894591) q[0];
rz(0.16667357) q[1];
sx q[1];
rz(-2.1107626) q[1];
sx q[1];
rz(-1.4717685) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6643065) q[0];
sx q[0];
rz(-2.1843231) q[0];
sx q[0];
rz(1.0856347) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6138863) q[2];
sx q[2];
rz(-2.913007) q[2];
sx q[2];
rz(3.0800386) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1907726) q[1];
sx q[1];
rz(-2.4476978) q[1];
sx q[1];
rz(-0.67339749) q[1];
rz(-pi) q[2];
rz(2.3237866) q[3];
sx q[3];
rz(-1.3997625) q[3];
sx q[3];
rz(0.59137646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9071831) q[2];
sx q[2];
rz(-3.0295591) q[2];
sx q[2];
rz(0.33667931) q[2];
rz(1.178721) q[3];
sx q[3];
rz(-1.4320635) q[3];
sx q[3];
rz(0.36693507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97230279) q[0];
sx q[0];
rz(-2.4827835) q[0];
sx q[0];
rz(1.4038053) q[0];
rz(2.9804969) q[1];
sx q[1];
rz(-2.1889071) q[1];
sx q[1];
rz(2.5162627) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0705497) q[0];
sx q[0];
rz(-0.82938507) q[0];
sx q[0];
rz(-2.2706288) q[0];
x q[1];
rz(-0.95809816) q[2];
sx q[2];
rz(-2.4110458) q[2];
sx q[2];
rz(2.036866) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4703964) q[1];
sx q[1];
rz(-2.0866864) q[1];
sx q[1];
rz(3.0615648) q[1];
rz(-pi) q[2];
rz(-0.13631374) q[3];
sx q[3];
rz(-1.9637311) q[3];
sx q[3];
rz(-0.64717347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8620341) q[2];
sx q[2];
rz(-1.4906078) q[2];
sx q[2];
rz(-3.1107235) q[2];
rz(-0.75012642) q[3];
sx q[3];
rz(-2.153219) q[3];
sx q[3];
rz(-1.8376902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52791643) q[0];
sx q[0];
rz(-0.91687098) q[0];
sx q[0];
rz(-3.0314639) q[0];
rz(0.88422173) q[1];
sx q[1];
rz(-1.5343816) q[1];
sx q[1];
rz(2.0933188) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3582402) q[0];
sx q[0];
rz(-1.8723346) q[0];
sx q[0];
rz(1.3025955) q[0];
rz(-3.0763392) q[2];
sx q[2];
rz(-1.835923) q[2];
sx q[2];
rz(1.2733151) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5267262) q[1];
sx q[1];
rz(-1.6294565) q[1];
sx q[1];
rz(-0.23414302) q[1];
rz(0.061125867) q[3];
sx q[3];
rz(-1.3406274) q[3];
sx q[3];
rz(-1.597061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5732164) q[2];
sx q[2];
rz(-1.6171075) q[2];
sx q[2];
rz(-3.1061843) q[2];
rz(-3.0951485) q[3];
sx q[3];
rz(-1.2848264) q[3];
sx q[3];
rz(0.057028381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7923918) q[0];
sx q[0];
rz(-0.89850315) q[0];
sx q[0];
rz(0.13353702) q[0];
rz(-1.7298493) q[1];
sx q[1];
rz(-2.0202507) q[1];
sx q[1];
rz(-2.0917361) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2757169) q[0];
sx q[0];
rz(-1.9483074) q[0];
sx q[0];
rz(1.6927682) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0303512) q[2];
sx q[2];
rz(-1.5644511) q[2];
sx q[2];
rz(-2.480577) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.12956086) q[1];
sx q[1];
rz(-0.82967463) q[1];
sx q[1];
rz(3.0864703) q[1];
rz(-0.61232845) q[3];
sx q[3];
rz(-0.48053743) q[3];
sx q[3];
rz(-0.2162741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.92020804) q[2];
sx q[2];
rz(-0.68156663) q[2];
sx q[2];
rz(1.7690915) q[2];
rz(-1.3051322) q[3];
sx q[3];
rz(-1.2499481) q[3];
sx q[3];
rz(0.60177747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9735499) q[0];
sx q[0];
rz(-1.9815227) q[0];
sx q[0];
rz(2.2450182) q[0];
rz(2.1606481) q[1];
sx q[1];
rz(-0.70471057) q[1];
sx q[1];
rz(0.14095813) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25280639) q[0];
sx q[0];
rz(-1.8060291) q[0];
sx q[0];
rz(2.1034408) q[0];
rz(-3.1257764) q[2];
sx q[2];
rz(-2.2770303) q[2];
sx q[2];
rz(3.0045403) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2796717) q[1];
sx q[1];
rz(-2.1057317) q[1];
sx q[1];
rz(1.5831854) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3178918) q[3];
sx q[3];
rz(-2.0440954) q[3];
sx q[3];
rz(-0.85609037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7732546) q[2];
sx q[2];
rz(-2.1835623) q[2];
sx q[2];
rz(-0.21993302) q[2];
rz(-0.95647612) q[3];
sx q[3];
rz(-0.82436162) q[3];
sx q[3];
rz(-2.1970356) q[3];
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
rz(-pi) q[0];
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
rz(-1.3750951) q[0];
sx q[0];
rz(-1.5909593) q[0];
sx q[0];
rz(2.5084556) q[0];
rz(-1.7924869) q[1];
sx q[1];
rz(-2.2213885) q[1];
sx q[1];
rz(2.3470168) q[1];
rz(1.9021992) q[2];
sx q[2];
rz(-2.080658) q[2];
sx q[2];
rz(-1.5823564) q[2];
rz(2.9969146) q[3];
sx q[3];
rz(-2.288648) q[3];
sx q[3];
rz(-0.6821117) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
