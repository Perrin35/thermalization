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
rz(0.14468004) q[0];
rz(3.559685) q[1];
sx q[1];
rz(3.2453645) q[1];
sx q[1];
rz(11.054463) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1967752) q[0];
sx q[0];
rz(-2.631105) q[0];
sx q[0];
rz(-1.8148242) q[0];
x q[1];
rz(1.51654) q[2];
sx q[2];
rz(-2.3035604) q[2];
sx q[2];
rz(2.4861479) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6408893) q[1];
sx q[1];
rz(-1.2391866) q[1];
sx q[1];
rz(1.9808116) q[1];
x q[2];
rz(-2.6295794) q[3];
sx q[3];
rz(-1.1992992) q[3];
sx q[3];
rz(-0.12646261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8251553) q[2];
sx q[2];
rz(-1.8708159) q[2];
sx q[2];
rz(-0.6553418) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69219387) q[0];
sx q[0];
rz(-1.3757502) q[0];
sx q[0];
rz(-0.50814116) q[0];
rz(-1.7610158) q[1];
sx q[1];
rz(-0.5813798) q[1];
sx q[1];
rz(-1.9320206) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66099453) q[0];
sx q[0];
rz(-0.61836243) q[0];
sx q[0];
rz(1.8905276) q[0];
x q[1];
rz(0.31866535) q[2];
sx q[2];
rz(-3.0851611) q[2];
sx q[2];
rz(2.3669764) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0712378) q[1];
sx q[1];
rz(-2.6609328) q[1];
sx q[1];
rz(2.8621469) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0367514) q[3];
sx q[3];
rz(-2.4656714) q[3];
sx q[3];
rz(0.40646857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9380583) q[2];
sx q[2];
rz(-1.7455696) q[2];
sx q[2];
rz(-2.5246942) q[2];
rz(1.5361702) q[3];
sx q[3];
rz(-2.0974396) q[3];
sx q[3];
rz(-0.49083403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1598375) q[0];
sx q[0];
rz(-1.497739) q[0];
sx q[0];
rz(2.4165261) q[0];
rz(-1.2695351) q[1];
sx q[1];
rz(-2.4639362) q[1];
sx q[1];
rz(0.95917541) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77391439) q[0];
sx q[0];
rz(-1.1836021) q[0];
sx q[0];
rz(-1.0605124) q[0];
rz(0.41590642) q[2];
sx q[2];
rz(-1.9291483) q[2];
sx q[2];
rz(-0.22474536) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.37466418) q[1];
sx q[1];
rz(-2.9247724) q[1];
sx q[1];
rz(1.0105074) q[1];
x q[2];
rz(2.3183075) q[3];
sx q[3];
rz(-0.24634493) q[3];
sx q[3];
rz(2.5096517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3265257) q[2];
sx q[2];
rz(-1.2667789) q[2];
sx q[2];
rz(-0.25700021) q[2];
rz(2.0604996) q[3];
sx q[3];
rz(-0.5210146) q[3];
sx q[3];
rz(1.5787554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0839888) q[0];
sx q[0];
rz(-2.8631518) q[0];
sx q[0];
rz(1.9637015) q[0];
rz(0.15011694) q[1];
sx q[1];
rz(-2.6094486) q[1];
sx q[1];
rz(1.6252801) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1105791) q[0];
sx q[0];
rz(-0.4415126) q[0];
sx q[0];
rz(-2.2792321) q[0];
rz(-2.2994735) q[2];
sx q[2];
rz(-1.8531905) q[2];
sx q[2];
rz(-0.51674622) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1399556) q[1];
sx q[1];
rz(-1.4712584) q[1];
sx q[1];
rz(2.8147354) q[1];
rz(-pi) q[2];
rz(-0.64370207) q[3];
sx q[3];
rz(-1.4182404) q[3];
sx q[3];
rz(1.8495454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84245044) q[2];
sx q[2];
rz(-0.56988847) q[2];
sx q[2];
rz(2.2184856) q[2];
rz(-2.182377) q[3];
sx q[3];
rz(-2.1289289) q[3];
sx q[3];
rz(0.21624163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8206896) q[0];
sx q[0];
rz(-0.90383363) q[0];
sx q[0];
rz(0.88291105) q[0];
rz(-1.2954905) q[1];
sx q[1];
rz(-2.3675282) q[1];
sx q[1];
rz(-0.49497089) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9920693) q[0];
sx q[0];
rz(-1.8599293) q[0];
sx q[0];
rz(-1.2393214) q[0];
rz(-0.60762824) q[2];
sx q[2];
rz(-2.3403185) q[2];
sx q[2];
rz(-0.60333383) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5589757) q[1];
sx q[1];
rz(-2.4673273) q[1];
sx q[1];
rz(0.14087454) q[1];
rz(-pi) q[2];
rz(-0.35843973) q[3];
sx q[3];
rz(-1.0762941) q[3];
sx q[3];
rz(2.2089437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8098658) q[2];
sx q[2];
rz(-3.0220384) q[2];
sx q[2];
rz(-1.2681819) q[2];
rz(-3.0590893) q[3];
sx q[3];
rz(-1.637633) q[3];
sx q[3];
rz(2.4845607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6562011) q[0];
sx q[0];
rz(-2.9347561) q[0];
sx q[0];
rz(1.50151) q[0];
rz(-1.062475) q[1];
sx q[1];
rz(-1.1524009) q[1];
sx q[1];
rz(-1.1753488) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51477434) q[0];
sx q[0];
rz(-1.8747703) q[0];
sx q[0];
rz(-1.8016812) q[0];
x q[1];
rz(2.7926867) q[2];
sx q[2];
rz(-0.72427536) q[2];
sx q[2];
rz(1.5021715) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2986002) q[1];
sx q[1];
rz(-1.6879884) q[1];
sx q[1];
rz(2.8882364) q[1];
rz(-pi) q[2];
rz(-1.8683942) q[3];
sx q[3];
rz(-2.3087325) q[3];
sx q[3];
rz(2.0119865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8106653) q[2];
sx q[2];
rz(-1.8570447) q[2];
sx q[2];
rz(2.1539099) q[2];
rz(2.2641613) q[3];
sx q[3];
rz(-2.8743447) q[3];
sx q[3];
rz(-0.73006829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1838609) q[0];
sx q[0];
rz(-1.5393625) q[0];
sx q[0];
rz(-0.13885942) q[0];
rz(1.214437) q[1];
sx q[1];
rz(-1.6338394) q[1];
sx q[1];
rz(0.81659281) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3735376) q[0];
sx q[0];
rz(-1.5898503) q[0];
sx q[0];
rz(1.7956177) q[0];
x q[1];
rz(-0.76824378) q[2];
sx q[2];
rz(-1.7599987) q[2];
sx q[2];
rz(2.9598494) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5467756) q[1];
sx q[1];
rz(-0.71974126) q[1];
sx q[1];
rz(-1.627658) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2995629) q[3];
sx q[3];
rz(-2.2737164) q[3];
sx q[3];
rz(1.7563535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.82178086) q[2];
sx q[2];
rz(-2.767441) q[2];
sx q[2];
rz(-2.6569341) q[2];
rz(2.7759077) q[3];
sx q[3];
rz(-1.3644783) q[3];
sx q[3];
rz(-1.9827838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0115688) q[0];
sx q[0];
rz(-2.8304709) q[0];
sx q[0];
rz(0.097231641) q[0];
rz(-3.0367127) q[1];
sx q[1];
rz(-0.77969867) q[1];
sx q[1];
rz(1.8269151) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7044367) q[0];
sx q[0];
rz(-1.5771241) q[0];
sx q[0];
rz(0.0023715677) q[0];
rz(2.9179321) q[2];
sx q[2];
rz(-1.8065435) q[2];
sx q[2];
rz(-3.0652114) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4215678) q[1];
sx q[1];
rz(-2.0889258) q[1];
sx q[1];
rz(-1.6314477) q[1];
rz(-pi) q[2];
rz(-2.1076074) q[3];
sx q[3];
rz(-1.4958463) q[3];
sx q[3];
rz(2.1225464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9647727) q[2];
sx q[2];
rz(-1.3975881) q[2];
sx q[2];
rz(0.1532661) q[2];
rz(1.9249453) q[3];
sx q[3];
rz(-2.9235268) q[3];
sx q[3];
rz(0.6161859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1185054) q[0];
sx q[0];
rz(-0.0054792976) q[0];
sx q[0];
rz(-1.5047005) q[0];
rz(0.79832375) q[1];
sx q[1];
rz(-1.5232122) q[1];
sx q[1];
rz(0.33531478) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3753877) q[0];
sx q[0];
rz(-1.6339301) q[0];
sx q[0];
rz(1.8324018) q[0];
rz(-0.7503926) q[2];
sx q[2];
rz(-1.5529035) q[2];
sx q[2];
rz(0.6070348) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3533201) q[1];
sx q[1];
rz(-1.2306552) q[1];
sx q[1];
rz(-1.2659094) q[1];
rz(0.15729372) q[3];
sx q[3];
rz(-2.1424737) q[3];
sx q[3];
rz(1.8974822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0140784) q[2];
sx q[2];
rz(-1.9731584) q[2];
sx q[2];
rz(-0.43759313) q[2];
rz(2.0994999) q[3];
sx q[3];
rz(-1.7362678) q[3];
sx q[3];
rz(2.5991345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9893148) q[0];
sx q[0];
rz(-0.88459009) q[0];
sx q[0];
rz(0.44961318) q[0];
rz(0.97688976) q[1];
sx q[1];
rz(-1.6376817) q[1];
sx q[1];
rz(0.50122112) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6527443) q[0];
sx q[0];
rz(-0.34372675) q[0];
sx q[0];
rz(0.97596844) q[0];
rz(1.1193163) q[2];
sx q[2];
rz(-2.6767133) q[2];
sx q[2];
rz(-1.3550188) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7536555) q[1];
sx q[1];
rz(-0.5273312) q[1];
sx q[1];
rz(1.0868294) q[1];
x q[2];
rz(0.34256713) q[3];
sx q[3];
rz(-1.1293355) q[3];
sx q[3];
rz(1.6395237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6476743) q[0];
sx q[0];
rz(-1.7128581) q[0];
sx q[0];
rz(0.38059522) q[0];
rz(-2.4299798) q[1];
sx q[1];
rz(-1.8743534) q[1];
sx q[1];
rz(1.4030917) q[1];
rz(-3.0416476) q[2];
sx q[2];
rz(-2.6700085) q[2];
sx q[2];
rz(2.811583) q[2];
rz(1.8071411) q[3];
sx q[3];
rz(-1.9665039) q[3];
sx q[3];
rz(2.891091) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
