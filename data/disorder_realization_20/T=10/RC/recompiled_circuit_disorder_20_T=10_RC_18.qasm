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
rz(1.8600872) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6706657) q[0];
sx q[0];
rz(-2.7424194) q[0];
sx q[0];
rz(-2.1098718) q[0];
rz(-pi) q[1];
rz(-2.4490657) q[2];
sx q[2];
rz(-1.5480969) q[2];
sx q[2];
rz(1.9164617) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.68583381) q[1];
sx q[1];
rz(-1.403192) q[1];
sx q[1];
rz(1.140825) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.063135677) q[3];
sx q[3];
rz(-1.4423443) q[3];
sx q[3];
rz(1.1032784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.32221258) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(1.5167351) q[2];
rz(-0.24762282) q[3];
sx q[3];
rz(-1.6678436) q[3];
sx q[3];
rz(-0.65555278) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6973998) q[0];
sx q[0];
rz(-1.4780761) q[0];
sx q[0];
rz(2.9852988) q[0];
rz(-2.5022751) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(1.6527536) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3296559) q[0];
sx q[0];
rz(-0.57528472) q[0];
sx q[0];
rz(0.29559691) q[0];
x q[1];
rz(0.70398753) q[2];
sx q[2];
rz(-2.1805602) q[2];
sx q[2];
rz(-3.0733382) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9298047) q[1];
sx q[1];
rz(-1.025052) q[1];
sx q[1];
rz(0.8916698) q[1];
x q[2];
rz(0.21807166) q[3];
sx q[3];
rz(-3.0017188) q[3];
sx q[3];
rz(-1.1034031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8504101) q[2];
sx q[2];
rz(-0.098878421) q[2];
sx q[2];
rz(2.8448811) q[2];
rz(-0.3324278) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(-1.9077574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333703) q[0];
sx q[0];
rz(-1.8828266) q[0];
sx q[0];
rz(-2.6112774) q[0];
rz(-2.5976394) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(1.189032) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9246763) q[0];
sx q[0];
rz(-2.8042256) q[0];
sx q[0];
rz(2.6485373) q[0];
rz(-pi) q[1];
rz(0.15709917) q[2];
sx q[2];
rz(-0.25114775) q[2];
sx q[2];
rz(-2.7709393) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5570453) q[1];
sx q[1];
rz(-1.0439596) q[1];
sx q[1];
rz(-0.79750632) q[1];
rz(1.1071113) q[3];
sx q[3];
rz(-0.92671219) q[3];
sx q[3];
rz(2.7657699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90551463) q[2];
sx q[2];
rz(-1.6685852) q[2];
sx q[2];
rz(-1.8900324) q[2];
rz(0.45423147) q[3];
sx q[3];
rz(-2.7082704) q[3];
sx q[3];
rz(-2.7864578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3635451) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(1.3409412) q[0];
rz(-2.5355133) q[1];
sx q[1];
rz(-0.93901912) q[1];
sx q[1];
rz(1.1846503) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.143648) q[0];
sx q[0];
rz(-2.8638683) q[0];
sx q[0];
rz(-0.68302897) q[0];
x q[1];
rz(-2.0834288) q[2];
sx q[2];
rz(-1.1626273) q[2];
sx q[2];
rz(-2.4525016) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9222316) q[1];
sx q[1];
rz(-2.1582099) q[1];
sx q[1];
rz(-0.14031336) q[1];
rz(-pi) q[2];
rz(-1.9373478) q[3];
sx q[3];
rz(-1.1712211) q[3];
sx q[3];
rz(-1.2031581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.23345315) q[2];
sx q[2];
rz(-0.88580695) q[2];
sx q[2];
rz(-0.63956368) q[2];
rz(0.8574287) q[3];
sx q[3];
rz(-1.144751) q[3];
sx q[3];
rz(2.6517984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5109167) q[0];
sx q[0];
rz(-0.54172051) q[0];
sx q[0];
rz(2.6211991) q[0];
rz(0.10294542) q[1];
sx q[1];
rz(-2.091145) q[1];
sx q[1];
rz(-0.72881126) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9279014) q[0];
sx q[0];
rz(-1.3662845) q[0];
sx q[0];
rz(1.6708899) q[0];
rz(1.7545627) q[2];
sx q[2];
rz(-0.73256058) q[2];
sx q[2];
rz(1.5866605) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5506958) q[1];
sx q[1];
rz(-1.7753074) q[1];
sx q[1];
rz(-1.900161) q[1];
rz(0.95727386) q[3];
sx q[3];
rz(-2.0271218) q[3];
sx q[3];
rz(-0.43382713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0430498) q[2];
sx q[2];
rz(-2.7816732) q[2];
sx q[2];
rz(-0.70971242) q[2];
rz(2.5085311) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(0.13171296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0387886) q[0];
sx q[0];
rz(-0.10237256) q[0];
sx q[0];
rz(-2.2578755) q[0];
rz(-2.2209514) q[1];
sx q[1];
rz(-2.0795152) q[1];
sx q[1];
rz(-0.19764915) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55006856) q[0];
sx q[0];
rz(-2.1854489) q[0];
sx q[0];
rz(2.385925) q[0];
rz(-pi) q[1];
rz(-0.63981668) q[2];
sx q[2];
rz(-3.0090927) q[2];
sx q[2];
rz(-0.25623955) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.16285322) q[1];
sx q[1];
rz(-2.0254685) q[1];
sx q[1];
rz(1.6182369) q[1];
rz(-2.4297653) q[3];
sx q[3];
rz(-0.89755745) q[3];
sx q[3];
rz(2.2664546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1137696) q[2];
sx q[2];
rz(-1.343507) q[2];
sx q[2];
rz(-0.075604288) q[2];
rz(-1.6522225) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(0.82908019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41912115) q[0];
sx q[0];
rz(-0.88547456) q[0];
sx q[0];
rz(-0.042073123) q[0];
rz(1.9627337) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(-1.9721608) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6299316) q[0];
sx q[0];
rz(-0.17782623) q[0];
sx q[0];
rz(1.1849665) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.85497) q[2];
sx q[2];
rz(-1.1948164) q[2];
sx q[2];
rz(2.3665731) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.7536105) q[3];
sx q[3];
rz(-1.531732) q[3];
sx q[3];
rz(-0.29400533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3183257) q[2];
sx q[2];
rz(-2.0264758) q[2];
sx q[2];
rz(-3.1271093) q[2];
rz(-1.621834) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(-2.5884132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13977215) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(-0.56234223) q[0];
rz(1.4315804) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(-2.6838578) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.021488) q[0];
sx q[0];
rz(-2.0459667) q[0];
sx q[0];
rz(3.0475463) q[0];
x q[1];
rz(-2.2051815) q[2];
sx q[2];
rz(-0.80447703) q[2];
sx q[2];
rz(0.22125439) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8843958) q[1];
sx q[1];
rz(-2.307502) q[1];
sx q[1];
rz(3.1222621) q[1];
rz(-2.7982513) q[3];
sx q[3];
rz(-1.0700873) q[3];
sx q[3];
rz(0.2688558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3056425) q[2];
sx q[2];
rz(-2.2017411) q[2];
sx q[2];
rz(0.28253728) q[2];
rz(-0.95747581) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(2.6628475) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95865059) q[0];
sx q[0];
rz(-2.4562953) q[0];
sx q[0];
rz(1.7392993) q[0];
rz(2.4422586) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(1.8797849) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14784797) q[0];
sx q[0];
rz(-1.1435259) q[0];
sx q[0];
rz(-1.8673613) q[0];
rz(-pi) q[1];
rz(2.8724573) q[2];
sx q[2];
rz(-0.67101523) q[2];
sx q[2];
rz(3.105643) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3952179) q[1];
sx q[1];
rz(-0.74133855) q[1];
sx q[1];
rz(-0.37903255) q[1];
rz(-pi) q[2];
rz(0.34187596) q[3];
sx q[3];
rz(-0.55594) q[3];
sx q[3];
rz(-0.13560175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4385779) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(1.4578488) q[2];
rz(-0.66649377) q[3];
sx q[3];
rz(-1.4278744) q[3];
sx q[3];
rz(0.75171793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89501971) q[0];
sx q[0];
rz(-0.81037766) q[0];
sx q[0];
rz(1.1599468) q[0];
rz(-0.054140422) q[1];
sx q[1];
rz(-1.4777947) q[1];
sx q[1];
rz(-2.0711526) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1150072) q[0];
sx q[0];
rz(-0.75782776) q[0];
sx q[0];
rz(1.3269781) q[0];
rz(-pi) q[1];
rz(0.038254914) q[2];
sx q[2];
rz(-0.80637156) q[2];
sx q[2];
rz(0.56626608) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.479446) q[1];
sx q[1];
rz(-1.8188735) q[1];
sx q[1];
rz(2.3723888) q[1];
x q[2];
rz(-0.6108547) q[3];
sx q[3];
rz(-2.8419552) q[3];
sx q[3];
rz(2.6328997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76876172) q[2];
sx q[2];
rz(-2.5925345) q[2];
sx q[2];
rz(1.5268415) q[2];
rz(-2.5620143) q[3];
sx q[3];
rz(-1.3468346) q[3];
sx q[3];
rz(-0.65210623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51331818) q[0];
sx q[0];
rz(-1.5341298) q[0];
sx q[0];
rz(-2.4232724) q[0];
rz(1.0340446) q[1];
sx q[1];
rz(-2.4517192) q[1];
sx q[1];
rz(2.413961) q[1];
rz(1.5288011) q[2];
sx q[2];
rz(-0.73245807) q[2];
sx q[2];
rz(3.0277638) q[2];
rz(-2.0680239) q[3];
sx q[3];
rz(-2.0484925) q[3];
sx q[3];
rz(-1.9852553) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
