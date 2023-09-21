OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5946755) q[0];
sx q[0];
rz(-1.0008873) q[0];
sx q[0];
rz(2.9291908) q[0];
rz(-2.4266333) q[1];
sx q[1];
rz(-0.7874878) q[1];
sx q[1];
rz(-1.2815055) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0467487) q[0];
sx q[0];
rz(-1.2307414) q[0];
sx q[0];
rz(-2.9283471) q[0];
rz(-pi) q[1];
rz(-2.4490657) q[2];
sx q[2];
rz(-1.5480969) q[2];
sx q[2];
rz(1.9164617) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9077397) q[1];
sx q[1];
rz(-0.45957652) q[1];
sx q[1];
rz(-1.9563667) q[1];
rz(-pi) q[2];
rz(1.4420907) q[3];
sx q[3];
rz(-1.5081815) q[3];
sx q[3];
rz(-2.6821729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8193801) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(1.5167351) q[2];
rz(0.24762282) q[3];
sx q[3];
rz(-1.473749) q[3];
sx q[3];
rz(-0.65555278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6973998) q[0];
sx q[0];
rz(-1.4780761) q[0];
sx q[0];
rz(2.9852988) q[0];
rz(2.5022751) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(1.4888391) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6328837) q[0];
sx q[0];
rz(-1.4116305) q[0];
sx q[0];
rz(-0.55523086) q[0];
rz(-pi) q[1];
rz(-2.4376051) q[2];
sx q[2];
rz(-2.1805602) q[2];
sx q[2];
rz(0.068254452) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3540346) q[1];
sx q[1];
rz(-0.84317849) q[1];
sx q[1];
rz(0.80227279) q[1];
rz(-pi) q[2];
rz(-1.5403455) q[3];
sx q[3];
rz(-1.4342562) q[3];
sx q[3];
rz(2.2583435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8504101) q[2];
sx q[2];
rz(-0.098878421) q[2];
sx q[2];
rz(-2.8448811) q[2];
rz(-0.3324278) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(-1.9077574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5082224) q[0];
sx q[0];
rz(-1.8828266) q[0];
sx q[0];
rz(2.6112774) q[0];
rz(0.5439533) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(-1.9525607) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8231359) q[0];
sx q[0];
rz(-1.4134777) q[0];
sx q[0];
rz(-0.29968981) q[0];
rz(-1.5306773) q[2];
sx q[2];
rz(-1.3228068) q[2];
sx q[2];
rz(2.9330394) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.51057286) q[1];
sx q[1];
rz(-0.90386183) q[1];
sx q[1];
rz(-2.2651947) q[1];
rz(-pi) q[2];
rz(1.1071113) q[3];
sx q[3];
rz(-2.2148805) q[3];
sx q[3];
rz(-2.7657699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.236078) q[2];
sx q[2];
rz(-1.4730075) q[2];
sx q[2];
rz(1.8900324) q[2];
rz(-0.45423147) q[3];
sx q[3];
rz(-0.43332228) q[3];
sx q[3];
rz(0.35513487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77804756) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(-1.3409412) q[0];
rz(-0.60607934) q[1];
sx q[1];
rz(-0.93901912) q[1];
sx q[1];
rz(-1.1846503) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8458762) q[0];
sx q[0];
rz(-1.7850951) q[0];
sx q[0];
rz(-1.3927668) q[0];
rz(0.46063936) q[2];
sx q[2];
rz(-1.1038291) q[2];
sx q[2];
rz(2.4796782) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9222316) q[1];
sx q[1];
rz(-2.1582099) q[1];
sx q[1];
rz(3.0012793) q[1];
rz(-pi) q[2];
rz(-0.70373669) q[3];
sx q[3];
rz(-0.53547137) q[3];
sx q[3];
rz(-2.717201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9081395) q[2];
sx q[2];
rz(-2.2557857) q[2];
sx q[2];
rz(-0.63956368) q[2];
rz(2.284164) q[3];
sx q[3];
rz(-1.144751) q[3];
sx q[3];
rz(0.48979428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.63067591) q[0];
sx q[0];
rz(-0.54172051) q[0];
sx q[0];
rz(2.6211991) q[0];
rz(0.10294542) q[1];
sx q[1];
rz(-2.091145) q[1];
sx q[1];
rz(2.4127814) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6729352) q[0];
sx q[0];
rz(-0.22738439) q[0];
sx q[0];
rz(2.6926281) q[0];
x q[1];
rz(2.2949218) q[2];
sx q[2];
rz(-1.6933105) q[2];
sx q[2];
rz(0.12144897) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5908969) q[1];
sx q[1];
rz(-1.7753074) q[1];
sx q[1];
rz(-1.900161) q[1];
rz(-pi) q[2];
rz(2.1843188) q[3];
sx q[3];
rz(-1.1144708) q[3];
sx q[3];
rz(-0.43382713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0430498) q[2];
sx q[2];
rz(-0.35991943) q[2];
sx q[2];
rz(2.4318802) q[2];
rz(-0.63306159) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(-3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1028041) q[0];
sx q[0];
rz(-3.0392201) q[0];
sx q[0];
rz(-0.88371712) q[0];
rz(0.9206413) q[1];
sx q[1];
rz(-1.0620774) q[1];
sx q[1];
rz(0.19764915) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6186002) q[0];
sx q[0];
rz(-0.97609659) q[0];
sx q[0];
rz(-2.3408875) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4913885) q[2];
sx q[2];
rz(-1.6769771) q[2];
sx q[2];
rz(2.2413145) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.16285322) q[1];
sx q[1];
rz(-1.1161242) q[1];
sx q[1];
rz(1.5233558) q[1];
x q[2];
rz(-0.75943767) q[3];
sx q[3];
rz(-2.1067838) q[3];
sx q[3];
rz(-2.939455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1137696) q[2];
sx q[2];
rz(-1.343507) q[2];
sx q[2];
rz(0.075604288) q[2];
rz(-1.6522225) q[3];
sx q[3];
rz(-2.7421156) q[3];
sx q[3];
rz(-0.82908019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.7224715) q[0];
sx q[0];
rz(-0.88547456) q[0];
sx q[0];
rz(0.042073123) q[0];
rz(-1.178859) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(1.1694318) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1202576) q[0];
sx q[0];
rz(-1.406167) q[0];
sx q[0];
rz(-0.067532587) q[0];
rz(-pi) q[1];
rz(-0.61755007) q[2];
sx q[2];
rz(-0.46717656) q[2];
sx q[2];
rz(0.10305603) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6199854) q[1];
sx q[1];
rz(-2.1072525) q[1];
sx q[1];
rz(1.9636088) q[1];
rz(-pi) q[2];
rz(-1.6129937) q[3];
sx q[3];
rz(-1.9584667) q[3];
sx q[3];
rz(-1.8488415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3183257) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(0.014483359) q[2];
rz(1.5197586) q[3];
sx q[3];
rz(-1.2967916) q[3];
sx q[3];
rz(2.5884132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0018205) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(-0.56234223) q[0];
rz(1.4315804) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(-2.6838578) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4938175) q[0];
sx q[0];
rz(-1.4871948) q[0];
sx q[0];
rz(-2.0477707) q[0];
rz(-2.2675603) q[2];
sx q[2];
rz(-1.1296141) q[2];
sx q[2];
rz(2.2639084) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8843958) q[1];
sx q[1];
rz(-0.83409062) q[1];
sx q[1];
rz(3.1222621) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1222955) q[3];
sx q[3];
rz(-2.5428452) q[3];
sx q[3];
rz(-2.7703354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83595014) q[2];
sx q[2];
rz(-2.2017411) q[2];
sx q[2];
rz(0.28253728) q[2];
rz(-2.1841168) q[3];
sx q[3];
rz(-1.3922393) q[3];
sx q[3];
rz(2.6628475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1829421) q[0];
sx q[0];
rz(-2.4562953) q[0];
sx q[0];
rz(1.4022934) q[0];
rz(-0.69933403) q[1];
sx q[1];
rz(-1.8383086) q[1];
sx q[1];
rz(-1.8797849) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48755074) q[0];
sx q[0];
rz(-2.6267509) q[0];
sx q[0];
rz(2.5709855) q[0];
x q[1];
rz(0.26913531) q[2];
sx q[2];
rz(-0.67101523) q[2];
sx q[2];
rz(0.035949635) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3952179) q[1];
sx q[1];
rz(-0.74133855) q[1];
sx q[1];
rz(0.37903255) q[1];
rz(-0.34187596) q[3];
sx q[3];
rz(-0.55594) q[3];
sx q[3];
rz(0.13560175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4385779) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(1.4578488) q[2];
rz(-2.4750989) q[3];
sx q[3];
rz(-1.7137182) q[3];
sx q[3];
rz(-2.3898747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2465729) q[0];
sx q[0];
rz(-2.331215) q[0];
sx q[0];
rz(1.9816459) q[0];
rz(-0.054140422) q[1];
sx q[1];
rz(-1.663798) q[1];
sx q[1];
rz(-1.0704401) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1150072) q[0];
sx q[0];
rz(-2.3837649) q[0];
sx q[0];
rz(-1.3269781) q[0];
rz(-pi) q[1];
rz(-2.3355867) q[2];
sx q[2];
rz(-1.5984048) q[2];
sx q[2];
rz(1.031014) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9995607) q[1];
sx q[1];
rz(-2.3107717) q[1];
sx q[1];
rz(1.9097411) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6108547) q[3];
sx q[3];
rz(-0.29963747) q[3];
sx q[3];
rz(-2.6328997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.76876172) q[2];
sx q[2];
rz(-0.54905811) q[2];
sx q[2];
rz(-1.5268415) q[2];
rz(2.5620143) q[3];
sx q[3];
rz(-1.3468346) q[3];
sx q[3];
rz(-2.4894864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51331818) q[0];
sx q[0];
rz(-1.6074629) q[0];
sx q[0];
rz(0.71832023) q[0];
rz(2.1075481) q[1];
sx q[1];
rz(-0.68987344) q[1];
sx q[1];
rz(-0.72763163) q[1];
rz(0.037739567) q[2];
sx q[2];
rz(-0.83913091) q[2];
sx q[2];
rz(-0.057374949) q[2];
rz(2.3970849) q[3];
sx q[3];
rz(-0.67529882) q[3];
sx q[3];
rz(-2.8532367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];