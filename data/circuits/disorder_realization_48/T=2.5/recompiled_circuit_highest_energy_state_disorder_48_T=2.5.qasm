OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.78955221176147) q[0];
sx q[0];
rz(3.96652505000169) q[0];
sx q[0];
rz(12.0316028356473) q[0];
rz(0.984002470970154) q[1];
sx q[1];
rz(0.505529554682322) q[1];
sx q[1];
rz(11.3043917179029) q[1];
cx q[1],q[0];
rz(-0.136848628520966) q[0];
sx q[0];
rz(1.83335772355134) q[0];
sx q[0];
rz(11.6189851522367) q[0];
rz(2.1106538772583) q[2];
sx q[2];
rz(1.01540199120576) q[2];
sx q[2];
rz(9.12681630849048) q[2];
cx q[2],q[1];
rz(3.78851127624512) q[1];
sx q[1];
rz(3.69134673674638) q[1];
sx q[1];
rz(7.95324454306766) q[1];
rz(2.03766536712646) q[3];
sx q[3];
rz(6.71383133729035) q[3];
sx q[3];
rz(10.3854120731275) q[3];
cx q[3],q[2];
rz(-3.39062476158142) q[2];
sx q[2];
rz(6.79555645783479) q[2];
sx q[2];
rz(7.94174883364841) q[2];
rz(0.28462165594101) q[3];
sx q[3];
rz(3.76483079989488) q[3];
sx q[3];
rz(11.5021917581479) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.86772179603577) q[0];
sx q[0];
rz(5.55645743210847) q[0];
sx q[0];
rz(12.5253629445951) q[0];
rz(-4.34925174713135) q[1];
sx q[1];
rz(2.91378465493257) q[1];
sx q[1];
rz(14.0270437955777) q[1];
cx q[1],q[0];
rz(1.49672937393188) q[0];
sx q[0];
rz(1.71641996701295) q[0];
sx q[0];
rz(12.5365643262784) q[0];
rz(1.00701141357422) q[2];
sx q[2];
rz(3.99263951380784) q[2];
sx q[2];
rz(10.6144959688108) q[2];
cx q[2],q[1];
rz(1.16653120517731) q[1];
sx q[1];
rz(4.67377880414064) q[1];
sx q[1];
rz(9.16818008422061) q[1];
rz(0.23433019220829) q[3];
sx q[3];
rz(2.43574920495088) q[3];
sx q[3];
rz(10.3076155543248) q[3];
cx q[3],q[2];
rz(0.665512084960938) q[2];
sx q[2];
rz(4.82186213334138) q[2];
sx q[2];
rz(7.67349908351108) q[2];
rz(0.0811854004859924) q[3];
sx q[3];
rz(5.6142636855417) q[3];
sx q[3];
rz(7.98887274264499) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.74710333347321) q[0];
sx q[0];
rz(3.15386268043006) q[0];
sx q[0];
rz(10.03485003709) q[0];
rz(1.31297814846039) q[1];
sx q[1];
rz(3.8372222503) q[1];
sx q[1];
rz(9.97539398669406) q[1];
cx q[1],q[0];
rz(0.0742155015468597) q[0];
sx q[0];
rz(4.06052157481248) q[0];
sx q[0];
rz(9.88205519913837) q[0];
rz(-0.87564891576767) q[2];
sx q[2];
rz(2.73104289372499) q[2];
sx q[2];
rz(10.9310747146527) q[2];
cx q[2],q[1];
rz(-1.0119343996048) q[1];
sx q[1];
rz(4.90062704880769) q[1];
sx q[1];
rz(9.7292168199937) q[1];
rz(0.154626205563545) q[3];
sx q[3];
rz(3.89274397690827) q[3];
sx q[3];
rz(8.99831963180705) q[3];
cx q[3],q[2];
rz(2.66374707221985) q[2];
sx q[2];
rz(3.3140258957916) q[2];
sx q[2];
rz(8.35839209555789) q[2];
rz(1.15866982936859) q[3];
sx q[3];
rz(4.52466955979402) q[3];
sx q[3];
rz(9.38273322432443) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.310241937637329) q[0];
sx q[0];
rz(3.75475183327729) q[0];
sx q[0];
rz(10.3418431639592) q[0];
rz(-0.655455827713013) q[1];
sx q[1];
rz(2.23234453995759) q[1];
sx q[1];
rz(9.81436741947337) q[1];
cx q[1],q[0];
rz(4.02149343490601) q[0];
sx q[0];
rz(4.52609828312928) q[0];
sx q[0];
rz(7.93710241316959) q[0];
rz(0.685048282146454) q[2];
sx q[2];
rz(1.11944022973115) q[2];
sx q[2];
rz(11.1719534158628) q[2];
cx q[2],q[1];
rz(-1.12341892719269) q[1];
sx q[1];
rz(4.68319991429383) q[1];
sx q[1];
rz(9.82610285877391) q[1];
rz(-0.728581368923187) q[3];
sx q[3];
rz(5.54401794274385) q[3];
sx q[3];
rz(12.2576019525449) q[3];
cx q[3],q[2];
rz(0.831708252429962) q[2];
sx q[2];
rz(4.08412322600419) q[2];
sx q[2];
rz(8.07246813773319) q[2];
rz(0.748121976852417) q[3];
sx q[3];
rz(4.44564476807649) q[3];
sx q[3];
rz(10.9514312505643) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.303695023059845) q[0];
sx q[0];
rz(5.71336999733979) q[0];
sx q[0];
rz(10.7662178039472) q[0];
rz(0.981264233589172) q[1];
sx q[1];
rz(5.00186065037782) q[1];
sx q[1];
rz(11.856415963165) q[1];
cx q[1],q[0];
rz(1.9729346036911) q[0];
sx q[0];
rz(6.6095394213968) q[0];
sx q[0];
rz(7.9136658668439) q[0];
rz(-2.56548619270325) q[2];
sx q[2];
rz(4.9135903437906) q[2];
sx q[2];
rz(5.51516959666415) q[2];
cx q[2],q[1];
rz(-0.0404742136597633) q[1];
sx q[1];
rz(4.92484894593293) q[1];
sx q[1];
rz(11.1697519779126) q[1];
rz(0.0949166417121887) q[3];
sx q[3];
rz(2.16649845440919) q[3];
sx q[3];
rz(7.49242589472934) q[3];
cx q[3],q[2];
rz(2.18334102630615) q[2];
sx q[2];
rz(5.52285090287263) q[2];
sx q[2];
rz(11.6596579313199) q[2];
rz(-0.182361543178558) q[3];
sx q[3];
rz(1.44011083443696) q[3];
sx q[3];
rz(10.2832072138707) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.3458503484726) q[0];
sx q[0];
rz(4.09632817109162) q[0];
sx q[0];
rz(9.22305042146846) q[0];
rz(-1.78511035442352) q[1];
sx q[1];
rz(2.4701665361696) q[1];
sx q[1];
rz(10.5759304523389) q[1];
cx q[1],q[0];
rz(1.19304525852203) q[0];
sx q[0];
rz(6.7616664489084) q[0];
sx q[0];
rz(8.95834822057887) q[0];
rz(0.949988842010498) q[2];
sx q[2];
rz(0.706389578180858) q[2];
sx q[2];
rz(9.71155769228145) q[2];
cx q[2],q[1];
rz(-1.50733017921448) q[1];
sx q[1];
rz(4.31019488175447) q[1];
sx q[1];
rz(10.3711970209996) q[1];
rz(-0.227898731827736) q[3];
sx q[3];
rz(1.07124749024446) q[3];
sx q[3];
rz(7.98528740405246) q[3];
cx q[3],q[2];
rz(-0.987289369106293) q[2];
sx q[2];
rz(3.72809049685533) q[2];
sx q[2];
rz(9.07732260822459) q[2];
rz(-0.821849822998047) q[3];
sx q[3];
rz(4.91785863240296) q[3];
sx q[3];
rz(10.9486679792325) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.80079412460327) q[0];
sx q[0];
rz(3.78844681580598) q[0];
sx q[0];
rz(8.30156061648532) q[0];
rz(-2.71282911300659) q[1];
sx q[1];
rz(2.25184968312318) q[1];
sx q[1];
rz(9.57373779117271) q[1];
cx q[1],q[0];
rz(0.652294158935547) q[0];
sx q[0];
rz(0.517443807917186) q[0];
sx q[0];
rz(9.37992876245781) q[0];
rz(1.27477502822876) q[2];
sx q[2];
rz(1.37755766709382) q[2];
sx q[2];
rz(12.5018744230191) q[2];
cx q[2],q[1];
rz(-0.808544456958771) q[1];
sx q[1];
rz(2.36969486077363) q[1];
sx q[1];
rz(11.2082307100217) q[1];
rz(3.22288370132446) q[3];
sx q[3];
rz(3.54856920440728) q[3];
sx q[3];
rz(11.0371821880262) q[3];
cx q[3],q[2];
rz(4.15356969833374) q[2];
sx q[2];
rz(3.49954703648622) q[2];
sx q[2];
rz(6.60467193125888) q[2];
rz(1.11644113063812) q[3];
sx q[3];
rz(4.00047674973542) q[3];
sx q[3];
rz(10.8387455701749) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.15349340438843) q[0];
sx q[0];
rz(4.93859282334382) q[0];
sx q[0];
rz(9.40186653322681) q[0];
rz(-1.54943907260895) q[1];
sx q[1];
rz(1.04330817063386) q[1];
sx q[1];
rz(14.4955544233243) q[1];
cx q[1],q[0];
rz(0.592719733715057) q[0];
sx q[0];
rz(3.75195333560044) q[0];
sx q[0];
rz(8.95786691307231) q[0];
rz(-0.731051504611969) q[2];
sx q[2];
rz(4.43800464470918) q[2];
sx q[2];
rz(7.662101840965) q[2];
cx q[2],q[1];
rz(-1.35009944438934) q[1];
sx q[1];
rz(1.1731357892328) q[1];
sx q[1];
rz(12.8953506708066) q[1];
rz(0.713494300842285) q[3];
sx q[3];
rz(4.78522697289521) q[3];
sx q[3];
rz(8.23897037505313) q[3];
cx q[3],q[2];
rz(2.21080780029297) q[2];
sx q[2];
rz(2.94269404013688) q[2];
sx q[2];
rz(11.5690073728482) q[2];
rz(-1.25842440128326) q[3];
sx q[3];
rz(4.33259549935395) q[3];
sx q[3];
rz(10.655158019058) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.157994300127029) q[0];
sx q[0];
rz(2.48574343522126) q[0];
sx q[0];
rz(9.89915609954997) q[0];
rz(-3.6932635307312) q[1];
sx q[1];
rz(5.51469055016572) q[1];
sx q[1];
rz(3.79912850855991) q[1];
cx q[1],q[0];
rz(-1.04845869541168) q[0];
sx q[0];
rz(-0.0462997833914081) q[0];
sx q[0];
rz(7.31725142001315) q[0];
rz(1.40086340904236) q[2];
sx q[2];
rz(4.31176987488801) q[2];
sx q[2];
rz(7.42102859019443) q[2];
cx q[2],q[1];
rz(4.46188402175903) q[1];
sx q[1];
rz(1.58709147770936) q[1];
sx q[1];
rz(9.05639473199054) q[1];
rz(-0.14454673230648) q[3];
sx q[3];
rz(4.81643477280671) q[3];
sx q[3];
rz(9.76988772153064) q[3];
cx q[3],q[2];
rz(-0.23003301024437) q[2];
sx q[2];
rz(2.73427891929681) q[2];
sx q[2];
rz(11.2550685167234) q[2];
rz(0.714105486869812) q[3];
sx q[3];
rz(5.00084617932374) q[3];
sx q[3];
rz(11.99653003215) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.00622009346261621) q[0];
sx q[0];
rz(5.33397355874116) q[0];
sx q[0];
rz(8.81648454665347) q[0];
rz(-0.817812025547028) q[1];
sx q[1];
rz(1.96559706528718) q[1];
sx q[1];
rz(10.2577111482541) q[1];
cx q[1],q[0];
rz(-0.352276295423508) q[0];
sx q[0];
rz(6.70763746102388) q[0];
sx q[0];
rz(10.8698590755384) q[0];
rz(-0.256621420383453) q[2];
sx q[2];
rz(4.40790429909761) q[2];
sx q[2];
rz(13.7862500905912) q[2];
cx q[2],q[1];
rz(2.2257833480835) q[1];
sx q[1];
rz(1.43810096581513) q[1];
sx q[1];
rz(7.27021834849521) q[1];
rz(0.0313519909977913) q[3];
sx q[3];
rz(5.30521074135835) q[3];
sx q[3];
rz(10.9415715694348) q[3];
cx q[3],q[2];
rz(4.32446908950806) q[2];
sx q[2];
rz(2.78878379066522) q[2];
sx q[2];
rz(5.69703361987277) q[2];
rz(2.23852491378784) q[3];
sx q[3];
rz(3.7640552838617) q[3];
sx q[3];
rz(12.4805254697721) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.863908171653748) q[0];
sx q[0];
rz(3.89049384196336) q[0];
sx q[0];
rz(10.1016231536786) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-2.82028388977051) q[1];
sx q[1];
rz(1.5317361672693) q[1];
sx q[1];
rz(14.3866820096891) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-5.10722208023071) q[2];
sx q[2];
rz(4.70511880715425) q[2];
sx q[2];
rz(10.6306925773542) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.76878559589386) q[3];
sx q[3];
rz(5.79901734192903) q[3];
sx q[3];
rz(10.3785224318425) q[3];
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
