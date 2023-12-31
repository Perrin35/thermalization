OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.5469172) q[0];
sx q[0];
rz(4.1424799) q[0];
sx q[0];
rz(9.6371798) q[0];
rz(-2.4266333) q[1];
sx q[1];
rz(-0.7874878) q[1];
sx q[1];
rz(-1.2815055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0467487) q[0];
sx q[0];
rz(-1.9108512) q[0];
sx q[0];
rz(-2.9283471) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6002866) q[2];
sx q[2];
rz(-2.2631096) q[2];
sx q[2];
rz(0.32683795) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68583381) q[1];
sx q[1];
rz(-1.403192) q[1];
sx q[1];
rz(-1.140825) q[1];
rz(-pi) q[2];
rz(-1.699502) q[3];
sx q[3];
rz(-1.5081815) q[3];
sx q[3];
rz(-2.6821729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.32221258) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(-1.6248576) q[2];
rz(-2.8939698) q[3];
sx q[3];
rz(-1.473749) q[3];
sx q[3];
rz(-0.65555278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6973998) q[0];
sx q[0];
rz(-1.6635165) q[0];
sx q[0];
rz(2.9852988) q[0];
rz(-2.5022751) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(1.6527536) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1600906) q[0];
sx q[0];
rz(-2.1182051) q[0];
sx q[0];
rz(-1.3840958) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.312617) q[2];
sx q[2];
rz(-1.0114705) q[2];
sx q[2];
rz(1.0499357) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78755806) q[1];
sx q[1];
rz(-2.2984142) q[1];
sx q[1];
rz(-2.3393199) q[1];
rz(1.5403455) q[3];
sx q[3];
rz(-1.4342562) q[3];
sx q[3];
rz(0.88324916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.8504101) q[2];
sx q[2];
rz(-3.0427142) q[2];
sx q[2];
rz(2.8448811) q[2];
rz(0.3324278) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(-1.2338352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5082224) q[0];
sx q[0];
rz(-1.2587661) q[0];
sx q[0];
rz(2.6112774) q[0];
rz(-0.5439533) q[1];
sx q[1];
rz(-2.0201611) q[1];
sx q[1];
rz(-1.9525607) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2169164) q[0];
sx q[0];
rz(-2.8042256) q[0];
sx q[0];
rz(-2.6485373) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5306773) q[2];
sx q[2];
rz(-1.8187858) q[2];
sx q[2];
rz(-0.2085533) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4424858) q[1];
sx q[1];
rz(-2.2190296) q[1];
sx q[1];
rz(-2.4590765) q[1];
rz(-pi) q[2];
rz(1.1071113) q[3];
sx q[3];
rz(-0.92671219) q[3];
sx q[3];
rz(-0.37582276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90551463) q[2];
sx q[2];
rz(-1.6685852) q[2];
sx q[2];
rz(1.8900324) q[2];
rz(-2.6873612) q[3];
sx q[3];
rz(-0.43332228) q[3];
sx q[3];
rz(-0.35513487) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3635451) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(1.8006515) q[0];
rz(-2.5355133) q[1];
sx q[1];
rz(-0.93901912) q[1];
sx q[1];
rz(-1.9569424) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8458762) q[0];
sx q[0];
rz(-1.7850951) q[0];
sx q[0];
rz(1.3927668) q[0];
rz(-pi) q[1];
rz(2.6809533) q[2];
sx q[2];
rz(-1.1038291) q[2];
sx q[2];
rz(-2.4796782) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.21936101) q[1];
sx q[1];
rz(-0.9833828) q[1];
sx q[1];
rz(-3.0012793) q[1];
x q[2];
rz(2.7167927) q[3];
sx q[3];
rz(-1.2343171) q[3];
sx q[3];
rz(2.6257023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.23345315) q[2];
sx q[2];
rz(-0.88580695) q[2];
sx q[2];
rz(-2.502029) q[2];
rz(-0.8574287) q[3];
sx q[3];
rz(-1.144751) q[3];
sx q[3];
rz(0.48979428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63067591) q[0];
sx q[0];
rz(-2.5998721) q[0];
sx q[0];
rz(-2.6211991) q[0];
rz(0.10294542) q[1];
sx q[1];
rz(-2.091145) q[1];
sx q[1];
rz(2.4127814) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9279014) q[0];
sx q[0];
rz(-1.7753082) q[0];
sx q[0];
rz(-1.4707028) q[0];
rz(-pi) q[1];
rz(-0.16291933) q[2];
sx q[2];
rz(-0.85328057) q[2];
sx q[2];
rz(-1.7999072) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5908969) q[1];
sx q[1];
rz(-1.3662852) q[1];
sx q[1];
rz(1.900161) q[1];
rz(2.2767931) q[3];
sx q[3];
rz(-2.3949361) q[3];
sx q[3];
rz(1.4454696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0430498) q[2];
sx q[2];
rz(-0.35991943) q[2];
sx q[2];
rz(2.4318802) q[2];
rz(0.63306159) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0387886) q[0];
sx q[0];
rz(-3.0392201) q[0];
sx q[0];
rz(-2.2578755) q[0];
rz(-0.9206413) q[1];
sx q[1];
rz(-1.0620774) q[1];
sx q[1];
rz(-0.19764915) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5915241) q[0];
sx q[0];
rz(-2.1854489) q[0];
sx q[0];
rz(-0.7556677) q[0];
rz(-pi) q[1];
rz(-0.63981668) q[2];
sx q[2];
rz(-3.0090927) q[2];
sx q[2];
rz(-0.25623955) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3870961) q[1];
sx q[1];
rz(-1.5281786) q[1];
sx q[1];
rz(0.45511647) q[1];
rz(2.382155) q[3];
sx q[3];
rz(-1.0348088) q[3];
sx q[3];
rz(2.939455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.027823042) q[2];
sx q[2];
rz(-1.343507) q[2];
sx q[2];
rz(-3.0659884) q[2];
rz(1.6522225) q[3];
sx q[3];
rz(-2.7421156) q[3];
sx q[3];
rz(0.82908019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7224715) q[0];
sx q[0];
rz(-0.88547456) q[0];
sx q[0];
rz(0.042073123) q[0];
rz(-1.9627337) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(1.9721608) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4394546) q[0];
sx q[0];
rz(-1.5041782) q[0];
sx q[0];
rz(-1.735795) q[0];
x q[1];
rz(-2.7514236) q[2];
sx q[2];
rz(-1.3069659) q[2];
sx q[2];
rz(-2.2389776) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9387503) q[1];
sx q[1];
rz(-0.65333594) q[1];
sx q[1];
rz(2.5696239) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0386482) q[3];
sx q[3];
rz(-0.38984459) q[3];
sx q[3];
rz(-1.9600705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3183257) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(-3.1271093) q[2];
rz(1.5197586) q[3];
sx q[3];
rz(-1.2967916) q[3];
sx q[3];
rz(2.5884132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0018205) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(0.56234223) q[0];
rz(1.4315804) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(0.45773488) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91677529) q[0];
sx q[0];
rz(-2.6579034) q[0];
sx q[0];
rz(-1.7513357) q[0];
rz(-pi) q[1];
rz(-2.5896795) q[2];
sx q[2];
rz(-2.1898824) q[2];
sx q[2];
rz(2.1053095) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.815005) q[1];
sx q[1];
rz(-1.5564789) q[1];
sx q[1];
rz(0.83399764) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1222955) q[3];
sx q[3];
rz(-0.5987474) q[3];
sx q[3];
rz(0.37125722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3056425) q[2];
sx q[2];
rz(-0.93985158) q[2];
sx q[2];
rz(-2.8590554) q[2];
rz(-2.1841168) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(0.47874513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95865059) q[0];
sx q[0];
rz(-0.68529737) q[0];
sx q[0];
rz(1.4022934) q[0];
rz(2.4422586) q[1];
sx q[1];
rz(-1.8383086) q[1];
sx q[1];
rz(1.2618077) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6540419) q[0];
sx q[0];
rz(-2.6267509) q[0];
sx q[0];
rz(-0.57060711) q[0];
rz(0.26913531) q[2];
sx q[2];
rz(-0.67101523) q[2];
sx q[2];
rz(0.035949635) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7463747) q[1];
sx q[1];
rz(-0.74133855) q[1];
sx q[1];
rz(-0.37903255) q[1];
rz(0.34187596) q[3];
sx q[3];
rz(-0.55594) q[3];
sx q[3];
rz(3.0059909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4385779) q[2];
sx q[2];
rz(-1.6377178) q[2];
sx q[2];
rz(-1.6837439) q[2];
rz(2.4750989) q[3];
sx q[3];
rz(-1.4278744) q[3];
sx q[3];
rz(-2.3898747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2465729) q[0];
sx q[0];
rz(-2.331215) q[0];
sx q[0];
rz(1.1599468) q[0];
rz(-3.0874522) q[1];
sx q[1];
rz(-1.4777947) q[1];
sx q[1];
rz(-1.0704401) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0265855) q[0];
sx q[0];
rz(-2.3837649) q[0];
sx q[0];
rz(-1.8146145) q[0];
x q[1];
rz(-1.530933) q[2];
sx q[2];
rz(-0.76518744) q[2];
sx q[2];
rz(2.6305692) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4813303) q[1];
sx q[1];
rz(-2.3412625) q[1];
sx q[1];
rz(-2.7923613) q[1];
x q[2];
rz(2.8937267) q[3];
sx q[3];
rz(-1.7409179) q[3];
sx q[3];
rz(-2.6691013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3728309) q[2];
sx q[2];
rz(-0.54905811) q[2];
sx q[2];
rz(1.5268415) q[2];
rz(0.57957831) q[3];
sx q[3];
rz(-1.7947581) q[3];
sx q[3];
rz(0.65210623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.51331818) q[0];
sx q[0];
rz(-1.5341298) q[0];
sx q[0];
rz(-2.4232724) q[0];
rz(-1.0340446) q[1];
sx q[1];
rz(-0.68987344) q[1];
sx q[1];
rz(-0.72763163) q[1];
rz(-0.037739567) q[2];
sx q[2];
rz(-2.3024617) q[2];
sx q[2];
rz(3.0842177) q[2];
rz(-2.3970849) q[3];
sx q[3];
rz(-2.4662938) q[3];
sx q[3];
rz(0.288356) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
