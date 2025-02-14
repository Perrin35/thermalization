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
rz(-0.55627745) q[0];
sx q[0];
rz(3.102432) q[0];
sx q[0];
rz(10.089212) q[0];
rz(-0.18222624) q[1];
sx q[1];
rz(4.8290904) q[1];
sx q[1];
rz(11.157142) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1202224) q[0];
sx q[0];
rz(-2.5787368) q[0];
sx q[0];
rz(-2.8261306) q[0];
rz(-pi) q[1];
rz(2.6970942) q[2];
sx q[2];
rz(-1.8730361) q[2];
sx q[2];
rz(-1.5465178) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.67068995) q[1];
sx q[1];
rz(-1.8876612) q[1];
sx q[1];
rz(2.5554727) q[1];
rz(-pi) q[2];
rz(-1.1143941) q[3];
sx q[3];
rz(-0.39908394) q[3];
sx q[3];
rz(1.2582859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0613784) q[2];
sx q[2];
rz(-0.55157101) q[2];
sx q[2];
rz(-0.75458327) q[2];
rz(-1.6825698) q[3];
sx q[3];
rz(-2.2857234) q[3];
sx q[3];
rz(1.4065019) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.215312) q[0];
sx q[0];
rz(-2.5975241) q[0];
sx q[0];
rz(1.1625483) q[0];
rz(-2.4303719) q[1];
sx q[1];
rz(-2.4237207) q[1];
sx q[1];
rz(0.1444764) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4220336) q[0];
sx q[0];
rz(-1.8633909) q[0];
sx q[0];
rz(-1.4841311) q[0];
x q[1];
rz(0.8501022) q[2];
sx q[2];
rz(-2.0726225) q[2];
sx q[2];
rz(0.36862954) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3857984) q[1];
sx q[1];
rz(-0.82634631) q[1];
sx q[1];
rz(2.6296294) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1816335) q[3];
sx q[3];
rz(-2.4414821) q[3];
sx q[3];
rz(2.7896529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43313906) q[2];
sx q[2];
rz(-1.3044367) q[2];
sx q[2];
rz(-0.30678314) q[2];
rz(-0.77543801) q[3];
sx q[3];
rz(-0.38475761) q[3];
sx q[3];
rz(-0.14498372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48258346) q[0];
sx q[0];
rz(-2.6982396) q[0];
sx q[0];
rz(2.3495667) q[0];
rz(1.5576942) q[1];
sx q[1];
rz(-1.3521103) q[1];
sx q[1];
rz(2.1477594) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4531012) q[0];
sx q[0];
rz(-1.2198041) q[0];
sx q[0];
rz(-0.61601244) q[0];
x q[1];
rz(-2.5517919) q[2];
sx q[2];
rz(-1.4695393) q[2];
sx q[2];
rz(2.4112005) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.12056226) q[1];
sx q[1];
rz(-2.0018199) q[1];
sx q[1];
rz(0.31224183) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1937856) q[3];
sx q[3];
rz(-0.48156958) q[3];
sx q[3];
rz(1.6653614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1899679) q[2];
sx q[2];
rz(-0.8364532) q[2];
sx q[2];
rz(2.9659029) q[2];
rz(-0.22819337) q[3];
sx q[3];
rz(-0.56932813) q[3];
sx q[3];
rz(2.6364251) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93155414) q[0];
sx q[0];
rz(-2.2966972) q[0];
sx q[0];
rz(-1.5616052) q[0];
rz(-2.2700894) q[1];
sx q[1];
rz(-1.5305287) q[1];
sx q[1];
rz(2.9759488) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79194389) q[0];
sx q[0];
rz(-1.5719746) q[0];
sx q[0];
rz(-3.1399957) q[0];
x q[1];
rz(1.5890938) q[2];
sx q[2];
rz(-0.4457655) q[2];
sx q[2];
rz(0.95761567) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9479182) q[1];
sx q[1];
rz(-0.79371023) q[1];
sx q[1];
rz(3.0315184) q[1];
rz(-0.58067643) q[3];
sx q[3];
rz(-1.8044953) q[3];
sx q[3];
rz(-0.068598821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6988301) q[2];
sx q[2];
rz(-2.5593968) q[2];
sx q[2];
rz(-0.77787918) q[2];
rz(0.88464087) q[3];
sx q[3];
rz(-1.9886465) q[3];
sx q[3];
rz(0.9790498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1409461) q[0];
sx q[0];
rz(-1.0362754) q[0];
sx q[0];
rz(-2.272814) q[0];
rz(-2.2390305) q[1];
sx q[1];
rz(-1.9156009) q[1];
sx q[1];
rz(2.5696519) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7392699) q[0];
sx q[0];
rz(-2.2796541) q[0];
sx q[0];
rz(2.2174382) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4221572) q[2];
sx q[2];
rz(-2.509382) q[2];
sx q[2];
rz(-0.84144634) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8780788) q[1];
sx q[1];
rz(-1.1357508) q[1];
sx q[1];
rz(0.33532354) q[1];
x q[2];
rz(1.0056313) q[3];
sx q[3];
rz(-0.4400357) q[3];
sx q[3];
rz(-2.7081851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.14704412) q[2];
sx q[2];
rz(-0.93783718) q[2];
sx q[2];
rz(-2.5082972) q[2];
rz(-2.5607732) q[3];
sx q[3];
rz(-1.3588901) q[3];
sx q[3];
rz(2.4766428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-3.0112515) q[0];
sx q[0];
rz(-2.7524152) q[0];
sx q[0];
rz(-1.3379541) q[0];
rz(-1.7300216) q[1];
sx q[1];
rz(-0.4069702) q[1];
sx q[1];
rz(2.3445047) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4914843) q[0];
sx q[0];
rz(-2.8239125) q[0];
sx q[0];
rz(-1.0450715) q[0];
rz(-0.077812151) q[2];
sx q[2];
rz(-2.1350265) q[2];
sx q[2];
rz(-2.5112453) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89576777) q[1];
sx q[1];
rz(-1.8386158) q[1];
sx q[1];
rz(-0.16606776) q[1];
x q[2];
rz(0.18081801) q[3];
sx q[3];
rz(-2.083484) q[3];
sx q[3];
rz(3.0231383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7756614) q[2];
sx q[2];
rz(-2.4727827) q[2];
sx q[2];
rz(-0.68525165) q[2];
rz(-2.8819486) q[3];
sx q[3];
rz(-0.97340596) q[3];
sx q[3];
rz(1.9297622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.892136) q[0];
sx q[0];
rz(-2.9849755) q[0];
sx q[0];
rz(2.6020965) q[0];
rz(-2.9501713) q[1];
sx q[1];
rz(-1.4861264) q[1];
sx q[1];
rz(-0.44245455) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4726329) q[0];
sx q[0];
rz(-1.5704535) q[0];
sx q[0];
rz(-1.5799205) q[0];
rz(-0.39789756) q[2];
sx q[2];
rz(-2.1566763) q[2];
sx q[2];
rz(0.90998024) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.48414184) q[1];
sx q[1];
rz(-1.9062348) q[1];
sx q[1];
rz(-0.92745933) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78667647) q[3];
sx q[3];
rz(-1.0287675) q[3];
sx q[3];
rz(-1.4876798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.29360867) q[2];
sx q[2];
rz(-0.30985761) q[2];
sx q[2];
rz(0.7027182) q[2];
rz(-2.8637049) q[3];
sx q[3];
rz(-2.0746456) q[3];
sx q[3];
rz(1.8028629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8948995) q[0];
sx q[0];
rz(-2.2791857) q[0];
sx q[0];
rz(-0.53584164) q[0];
rz(0.72227532) q[1];
sx q[1];
rz(-1.7012137) q[1];
sx q[1];
rz(-0.78295082) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.10038) q[0];
sx q[0];
rz(-1.6579275) q[0];
sx q[0];
rz(-1.210733) q[0];
rz(1.6269685) q[2];
sx q[2];
rz(-1.1740285) q[2];
sx q[2];
rz(-2.6068991) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47470505) q[1];
sx q[1];
rz(-1.0586481) q[1];
sx q[1];
rz(0.94774232) q[1];
rz(2.2244455) q[3];
sx q[3];
rz(-0.65509812) q[3];
sx q[3];
rz(0.411471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7213514) q[2];
sx q[2];
rz(-2.6768117) q[2];
sx q[2];
rz(-0.50344023) q[2];
rz(-2.8069046) q[3];
sx q[3];
rz(-1.3379593) q[3];
sx q[3];
rz(-0.23505841) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6675785) q[0];
sx q[0];
rz(-2.773371) q[0];
sx q[0];
rz(2.0817122) q[0];
rz(-0.31479752) q[1];
sx q[1];
rz(-1.7960725) q[1];
sx q[1];
rz(1.559929) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81875694) q[0];
sx q[0];
rz(-2.4222932) q[0];
sx q[0];
rz(2.8000205) q[0];
rz(3.1078075) q[2];
sx q[2];
rz(-1.3063161) q[2];
sx q[2];
rz(2.1839301) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.25645721) q[1];
sx q[1];
rz(-0.68084913) q[1];
sx q[1];
rz(1.225994) q[1];
x q[2];
rz(1.6208036) q[3];
sx q[3];
rz(-0.41566089) q[3];
sx q[3];
rz(2.5396944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18610893) q[2];
sx q[2];
rz(-1.1649705) q[2];
sx q[2];
rz(-2.1060409) q[2];
rz(-1.1459076) q[3];
sx q[3];
rz(-2.0290387) q[3];
sx q[3];
rz(-2.9884393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.938852) q[0];
sx q[0];
rz(-0.11883141) q[0];
sx q[0];
rz(-0.76549292) q[0];
rz(-2.1047523) q[1];
sx q[1];
rz(-1.2988657) q[1];
sx q[1];
rz(-0.11229215) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0033007) q[0];
sx q[0];
rz(-1.1213741) q[0];
sx q[0];
rz(1.9339191) q[0];
rz(-2.5300171) q[2];
sx q[2];
rz(-0.87885746) q[2];
sx q[2];
rz(2.126978) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0794509) q[1];
sx q[1];
rz(-1.7711258) q[1];
sx q[1];
rz(0.82605861) q[1];
rz(-pi) q[2];
rz(-0.78386098) q[3];
sx q[3];
rz(-1.5376159) q[3];
sx q[3];
rz(-2.9110094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3969193) q[2];
sx q[2];
rz(-1.2770709) q[2];
sx q[2];
rz(0.12167682) q[2];
rz(-2.7820898) q[3];
sx q[3];
rz(-0.72327852) q[3];
sx q[3];
rz(-0.017688964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.551238) q[0];
sx q[0];
rz(-1.4689162) q[0];
sx q[0];
rz(1.3015251) q[0];
rz(1.3941258) q[1];
sx q[1];
rz(-1.196741) q[1];
sx q[1];
rz(-1.7653042) q[1];
rz(2.0040705) q[2];
sx q[2];
rz(-2.1480297) q[2];
sx q[2];
rz(0.99110023) q[2];
rz(2.0064374) q[3];
sx q[3];
rz(-1.451685) q[3];
sx q[3];
rz(1.4901037) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
