OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.528429269790649) q[0];
sx q[0];
rz(2.08187213738496) q[0];
sx q[0];
rz(11.8353974580686) q[0];
rz(1.64147400856018) q[1];
sx q[1];
rz(5.24836221535737) q[1];
sx q[1];
rz(11.622826075546) q[1];
cx q[1],q[0];
rz(0.659817278385162) q[0];
sx q[0];
rz(3.21718379308517) q[0];
sx q[0];
rz(8.33548221587344) q[0];
rz(1.68204355239868) q[2];
sx q[2];
rz(5.6292304118448) q[2];
sx q[2];
rz(5.92047951220676) q[2];
cx q[2],q[1];
rz(-4.04884099960327) q[1];
sx q[1];
rz(5.75307169755036) q[1];
sx q[1];
rz(10.1445289611737) q[1];
rz(3.32898139953613) q[3];
sx q[3];
rz(1.87057796319062) q[3];
sx q[3];
rz(9.41633212137177) q[3];
cx q[3],q[2];
rz(-0.255081832408905) q[2];
sx q[2];
rz(1.76048973401124) q[2];
sx q[2];
rz(11.3156046628873) q[2];
rz(1.42617738246918) q[3];
sx q[3];
rz(4.05766061146791) q[3];
sx q[3];
rz(11.5864298105161) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.00868666172028) q[0];
sx q[0];
rz(5.21089473565156) q[0];
sx q[0];
rz(8.82583311795398) q[0];
rz(1.34097456932068) q[1];
sx q[1];
rz(4.09180882771546) q[1];
sx q[1];
rz(10.3911775112073) q[1];
cx q[1],q[0];
rz(0.371062844991684) q[0];
sx q[0];
rz(3.90096250374848) q[0];
sx q[0];
rz(11.6636142492215) q[0];
rz(2.19298005104065) q[2];
sx q[2];
rz(4.67195740540559) q[2];
sx q[2];
rz(7.68687567710086) q[2];
cx q[2],q[1];
rz(-0.84345531463623) q[1];
sx q[1];
rz(4.82624140580232) q[1];
sx q[1];
rz(12.6459862947385) q[1];
rz(1.1130485534668) q[3];
sx q[3];
rz(3.94509157736833) q[3];
sx q[3];
rz(10.4018271922986) q[3];
cx q[3],q[2];
rz(-6.1975269317627) q[2];
sx q[2];
rz(3.95674702723558) q[2];
sx q[2];
rz(15.406867003433) q[2];
rz(1.19312334060669) q[3];
sx q[3];
rz(4.69171515305574) q[3];
sx q[3];
rz(10.6111635923307) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.24626493453979) q[0];
sx q[0];
rz(4.7788623889261) q[0];
sx q[0];
rz(10.6122416019361) q[0];
rz(-1.90564119815826) q[1];
sx q[1];
rz(5.24583736260469) q[1];
sx q[1];
rz(10.7423640251081) q[1];
cx q[1],q[0];
rz(-1.00619578361511) q[0];
sx q[0];
rz(5.5580622275644) q[0];
sx q[0];
rz(11.3186484336774) q[0];
rz(-0.197443500161171) q[2];
sx q[2];
rz(1.78064456780488) q[2];
sx q[2];
rz(10.3393128871839) q[2];
cx q[2],q[1];
rz(4.61460018157959) q[1];
sx q[1];
rz(5.85044589837129) q[1];
sx q[1];
rz(11.2390191316526) q[1];
rz(1.68812012672424) q[3];
sx q[3];
rz(2.08527198632295) q[3];
sx q[3];
rz(9.16020647286578) q[3];
cx q[3],q[2];
rz(2.01117706298828) q[2];
sx q[2];
rz(1.39355638821656) q[2];
sx q[2];
rz(12.3597345113675) q[2];
rz(0.708042800426483) q[3];
sx q[3];
rz(2.93386241992051) q[3];
sx q[3];
rz(8.53195360898181) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.774650394916534) q[0];
sx q[0];
rz(2.92705244024331) q[0];
sx q[0];
rz(8.53855702876254) q[0];
rz(2.13185024261475) q[1];
sx q[1];
rz(5.37703243096406) q[1];
sx q[1];
rz(11.3399470806043) q[1];
cx q[1],q[0];
rz(3.28178644180298) q[0];
sx q[0];
rz(1.63700190384919) q[0];
sx q[0];
rz(8.69735071658298) q[0];
rz(4.25252389907837) q[2];
sx q[2];
rz(5.34585657914216) q[2];
sx q[2];
rz(11.299705362312) q[2];
cx q[2],q[1];
rz(4.28504323959351) q[1];
sx q[1];
rz(1.36100235779817) q[1];
sx q[1];
rz(9.62609378098651) q[1];
rz(-1.08441531658173) q[3];
sx q[3];
rz(4.3796547969156) q[3];
sx q[3];
rz(8.51625767945453) q[3];
cx q[3],q[2];
rz(1.64407885074615) q[2];
sx q[2];
rz(4.45546201069886) q[2];
sx q[2];
rz(11.5733740091245) q[2];
rz(-1.82890605926514) q[3];
sx q[3];
rz(2.11951807339723) q[3];
sx q[3];
rz(9.90851364134952) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.11836695671082) q[0];
sx q[0];
rz(2.82677304943139) q[0];
sx q[0];
rz(10.6107657909314) q[0];
rz(-1.42333078384399) q[1];
sx q[1];
rz(7.9344050010019) q[1];
sx q[1];
rz(6.86566255091831) q[1];
cx q[1],q[0];
rz(-1.54861485958099) q[0];
sx q[0];
rz(0.606152208643504) q[0];
sx q[0];
rz(6.15714905261203) q[0];
rz(-1.42180967330933) q[2];
sx q[2];
rz(2.2431230862909) q[2];
sx q[2];
rz(9.00178534387752) q[2];
cx q[2],q[1];
rz(-5.87699747085571) q[1];
sx q[1];
rz(4.69319990475709) q[1];
sx q[1];
rz(11.7100727319638) q[1];
rz(-0.0670276060700417) q[3];
sx q[3];
rz(5.11785581906373) q[3];
sx q[3];
rz(10.4321126699369) q[3];
cx q[3],q[2];
rz(-3.79645848274231) q[2];
sx q[2];
rz(4.74836161931092) q[2];
sx q[2];
rz(12.6480815172116) q[2];
rz(-0.47406667470932) q[3];
sx q[3];
rz(4.40523031552369) q[3];
sx q[3];
rz(11.0678174257199) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.89664697647095) q[0];
sx q[0];
rz(5.15521803696687) q[0];
sx q[0];
rz(9.41697303055927) q[0];
rz(1.40049779415131) q[1];
sx q[1];
rz(3.99565336306626) q[1];
sx q[1];
rz(10.5294241666715) q[1];
cx q[1],q[0];
rz(2.45451259613037) q[0];
sx q[0];
rz(5.36990896065766) q[0];
sx q[0];
rz(9.2183089017789) q[0];
rz(-2.21092653274536) q[2];
sx q[2];
rz(4.90862432320649) q[2];
sx q[2];
rz(6.18920061587497) q[2];
cx q[2],q[1];
rz(4.11723375320435) q[1];
sx q[1];
rz(4.01725325186784) q[1];
sx q[1];
rz(5.45167062281772) q[1];
rz(0.465521693229675) q[3];
sx q[3];
rz(4.86165431340272) q[3];
sx q[3];
rz(10.0309116601865) q[3];
cx q[3],q[2];
rz(-1.25099873542786) q[2];
sx q[2];
rz(3.8712678869539) q[2];
sx q[2];
rz(9.98000833987399) q[2];
rz(0.172720775008202) q[3];
sx q[3];
rz(4.96967867215211) q[3];
sx q[3];
rz(10.9730386495511) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.58844876289368) q[0];
sx q[0];
rz(3.62344169815118) q[0];
sx q[0];
rz(9.36301150395676) q[0];
rz(-0.242083668708801) q[1];
sx q[1];
rz(-0.375290719670705) q[1];
sx q[1];
rz(11.4544870614926) q[1];
cx q[1],q[0];
rz(2.04308485984802) q[0];
sx q[0];
rz(0.831391723948069) q[0];
sx q[0];
rz(8.84892699717685) q[0];
rz(1.55697023868561) q[2];
sx q[2];
rz(4.63478103478486) q[2];
sx q[2];
rz(8.11759016513034) q[2];
cx q[2],q[1];
rz(1.89167869091034) q[1];
sx q[1];
rz(0.877124460535594) q[1];
sx q[1];
rz(9.0616858959119) q[1];
rz(0.935636639595032) q[3];
sx q[3];
rz(6.68343177636201) q[3];
sx q[3];
rz(8.28981945513889) q[3];
cx q[3],q[2];
rz(3.95093655586243) q[2];
sx q[2];
rz(5.51884666283662) q[2];
sx q[2];
rz(8.61113355158969) q[2];
rz(-1.73711967468262) q[3];
sx q[3];
rz(3.37029591401155) q[3];
sx q[3];
rz(8.80935875176593) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.62548816204071) q[0];
sx q[0];
rz(4.49186423619325) q[0];
sx q[0];
rz(11.1409215688626) q[0];
rz(-1.52159929275513) q[1];
sx q[1];
rz(5.53125992615754) q[1];
sx q[1];
rz(6.92070314883395) q[1];
cx q[1],q[0];
rz(2.22560477256775) q[0];
sx q[0];
rz(3.99105009635026) q[0];
sx q[0];
rz(8.61808702944919) q[0];
rz(0.357122391462326) q[2];
sx q[2];
rz(4.01999953587586) q[2];
sx q[2];
rz(10.5784098863523) q[2];
cx q[2],q[1];
rz(-2.43791651725769) q[1];
sx q[1];
rz(5.29103937943513) q[1];
sx q[1];
rz(9.90603832005664) q[1];
rz(2.14150524139404) q[3];
sx q[3];
rz(4.82829669316346) q[3];
sx q[3];
rz(10.8259987592618) q[3];
cx q[3],q[2];
rz(-2.13115382194519) q[2];
sx q[2];
rz(3.84178742964799) q[2];
sx q[2];
rz(13.7025408506314) q[2];
rz(1.48539590835571) q[3];
sx q[3];
rz(0.524320038157054) q[3];
sx q[3];
rz(10.6343178510587) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.62563252449036) q[0];
sx q[0];
rz(3.92332551081712) q[0];
sx q[0];
rz(10.1009251832883) q[0];
rz(-2.31620836257935) q[1];
sx q[1];
rz(2.71765795548493) q[1];
sx q[1];
rz(11.5511924981992) q[1];
cx q[1],q[0];
rz(-0.434682011604309) q[0];
sx q[0];
rz(5.84594574769075) q[0];
sx q[0];
rz(7.54154524802371) q[0];
rz(0.153765425086021) q[2];
sx q[2];
rz(5.35948649247224) q[2];
sx q[2];
rz(3.67873332499667) q[2];
cx q[2],q[1];
rz(-3.20113396644592) q[1];
sx q[1];
rz(3.64812246163423) q[1];
sx q[1];
rz(10.2101028919141) q[1];
rz(1.22086989879608) q[3];
sx q[3];
rz(-0.797779885930471) q[3];
sx q[3];
rz(7.91607568263217) q[3];
cx q[3],q[2];
rz(0.721458077430725) q[2];
sx q[2];
rz(6.49942770798738) q[2];
sx q[2];
rz(7.25053641795322) q[2];
rz(-1.54452764987946) q[3];
sx q[3];
rz(4.45440021355683) q[3];
sx q[3];
rz(9.777650243036) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.618514597415924) q[0];
sx q[0];
rz(1.56275621254975) q[0];
sx q[0];
rz(9.48127591087624) q[0];
rz(-2.06911945343018) q[1];
sx q[1];
rz(5.01344767411286) q[1];
sx q[1];
rz(10.8294011115949) q[1];
cx q[1],q[0];
rz(-0.79849761724472) q[0];
sx q[0];
rz(4.58326760132844) q[0];
sx q[0];
rz(9.06139958500072) q[0];
rz(-1.19910967350006) q[2];
sx q[2];
rz(5.95823231537873) q[2];
sx q[2];
rz(4.77068946360751) q[2];
cx q[2],q[1];
rz(-0.235525876283646) q[1];
sx q[1];
rz(4.248159917193) q[1];
sx q[1];
rz(11.8105318307798) q[1];
rz(-0.00874906592071056) q[3];
sx q[3];
rz(7.70345083077485) q[3];
sx q[3];
rz(9.94262204169437) q[3];
cx q[3],q[2];
rz(0.29356062412262) q[2];
sx q[2];
rz(3.89548656542832) q[2];
sx q[2];
rz(9.73090088962718) q[2];
rz(-0.398115783929825) q[3];
sx q[3];
rz(4.80697027047212) q[3];
sx q[3];
rz(10.4490076064984) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.48306548595428) q[0];
sx q[0];
rz(5.12903145154054) q[0];
sx q[0];
rz(9.98214170931979) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-3.40008330345154) q[1];
sx q[1];
rz(4.46893409092958) q[1];
sx q[1];
rz(11.9981353044431) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(1.40100753307343) q[2];
sx q[2];
rz(3.66040715773637) q[2];
sx q[2];
rz(9.95202425717517) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.912241518497467) q[3];
sx q[3];
rz(4.65439441998536) q[3];
sx q[3];
rz(8.41971335410281) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];