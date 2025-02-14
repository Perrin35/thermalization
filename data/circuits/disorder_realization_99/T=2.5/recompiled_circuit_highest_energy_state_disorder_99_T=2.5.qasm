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
rz(1.9361629486084) q[0];
sx q[0];
rz(3.08714502503211) q[0];
sx q[0];
rz(8.55893907546207) q[0];
rz(7.88534498214722) q[1];
sx q[1];
rz(1.24422338803346) q[1];
sx q[1];
rz(3.19984529017612) q[1];
cx q[1],q[0];
rz(3.25235795974731) q[0];
sx q[0];
rz(4.54805460770661) q[0];
sx q[0];
rz(7.90155098437473) q[0];
rz(-0.136138260364532) q[2];
sx q[2];
rz(5.78364959557588) q[2];
sx q[2];
rz(10.9220494985501) q[2];
cx q[2],q[1];
rz(-0.213197514414787) q[1];
sx q[1];
rz(1.00939122040803) q[1];
sx q[1];
rz(11.9278824090879) q[1];
rz(-2.79440760612488) q[3];
sx q[3];
rz(4.40681126912171) q[3];
sx q[3];
rz(7.67680523394748) q[3];
cx q[3],q[2];
rz(1.31293475627899) q[2];
sx q[2];
rz(4.02223584254319) q[2];
sx q[2];
rz(10.3890746593396) q[2];
rz(0.0794845521450043) q[3];
sx q[3];
rz(4.98054495652253) q[3];
sx q[3];
rz(13.3375513315122) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-3.12817788124084) q[0];
sx q[0];
rz(5.93977394898469) q[0];
sx q[0];
rz(10.9651683330457) q[0];
rz(0.841142892837524) q[1];
sx q[1];
rz(4.54956653912599) q[1];
sx q[1];
rz(11.7089342832486) q[1];
cx q[1],q[0];
rz(1.67331373691559) q[0];
sx q[0];
rz(0.438600929575511) q[0];
sx q[0];
rz(9.4375957282181) q[0];
rz(-2.50560259819031) q[2];
sx q[2];
rz(5.10271409352357) q[2];
sx q[2];
rz(10.6019195079724) q[2];
cx q[2],q[1];
rz(-1.40743601322174) q[1];
sx q[1];
rz(3.8498382290178) q[1];
sx q[1];
rz(11.7042708158414) q[1];
rz(1.84990251064301) q[3];
sx q[3];
rz(1.66309038003022) q[3];
sx q[3];
rz(11.0219384193341) q[3];
cx q[3],q[2];
rz(0.886058270931244) q[2];
sx q[2];
rz(3.38730623026425) q[2];
sx q[2];
rz(8.21516082286044) q[2];
rz(1.38131928443909) q[3];
sx q[3];
rz(5.02229502995545) q[3];
sx q[3];
rz(5.69296524523898) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.48635387420654) q[0];
sx q[0];
rz(1.82625702221925) q[0];
sx q[0];
rz(8.09259340762302) q[0];
rz(1.47651290893555) q[1];
sx q[1];
rz(4.47246542771394) q[1];
sx q[1];
rz(6.90298864840671) q[1];
cx q[1],q[0];
rz(3.56981110572815) q[0];
sx q[0];
rz(3.22049505461986) q[0];
sx q[0];
rz(10.7959056854169) q[0];
rz(1.73930466175079) q[2];
sx q[2];
rz(4.96446815331513) q[2];
sx q[2];
rz(7.60767993926212) q[2];
cx q[2],q[1];
rz(0.425927489995956) q[1];
sx q[1];
rz(8.17112317879731) q[1];
sx q[1];
rz(9.34951984732553) q[1];
rz(-0.639846980571747) q[3];
sx q[3];
rz(4.92652037938172) q[3];
sx q[3];
rz(10.2242721080701) q[3];
cx q[3],q[2];
rz(-0.000981311197392642) q[2];
sx q[2];
rz(2.81622138817842) q[2];
sx q[2];
rz(11.7780163049619) q[2];
rz(5.00704050064087) q[3];
sx q[3];
rz(5.05593887169892) q[3];
sx q[3];
rz(10.2875911950986) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.15144228935242) q[0];
sx q[0];
rz(1.38624146779115) q[0];
sx q[0];
rz(10.4914899825971) q[0];
rz(1.78812956809998) q[1];
sx q[1];
rz(5.41624203522737) q[1];
sx q[1];
rz(8.26846120356723) q[1];
cx q[1],q[0];
rz(4.47177600860596) q[0];
sx q[0];
rz(4.69359138806398) q[0];
sx q[0];
rz(9.28324109911128) q[0];
rz(-0.040751788765192) q[2];
sx q[2];
rz(4.73857823212678) q[2];
sx q[2];
rz(14.5013942480008) q[2];
cx q[2],q[1];
rz(-5.11575746536255) q[1];
sx q[1];
rz(7.56515565712983) q[1];
sx q[1];
rz(9.81430054306194) q[1];
rz(0.906406044960022) q[3];
sx q[3];
rz(4.66023913224275) q[3];
sx q[3];
rz(8.46309570073291) q[3];
cx q[3],q[2];
rz(0.689169645309448) q[2];
sx q[2];
rz(3.8049974163347) q[2];
sx q[2];
rz(9.48360325246259) q[2];
rz(3.78092002868652) q[3];
sx q[3];
rz(5.06181195576722) q[3];
sx q[3];
rz(7.14053103923007) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.04621291160583) q[0];
sx q[0];
rz(7.69863525231416) q[0];
sx q[0];
rz(9.41428733839794) q[0];
rz(4.72016429901123) q[1];
sx q[1];
rz(3.83226427634294) q[1];
sx q[1];
rz(8.93541798590823) q[1];
cx q[1],q[0];
rz(-2.68835139274597) q[0];
sx q[0];
rz(5.63883152802522) q[0];
sx q[0];
rz(7.58924553393527) q[0];
rz(0.478859424591064) q[2];
sx q[2];
rz(3.89140293200547) q[2];
sx q[2];
rz(6.13031527995273) q[2];
cx q[2],q[1];
rz(-2.60766792297363) q[1];
sx q[1];
rz(8.20673909981782) q[1];
sx q[1];
rz(10.3006664276044) q[1];
rz(-1.23946607112885) q[3];
sx q[3];
rz(3.17412818049128) q[3];
sx q[3];
rz(14.6506247281949) q[3];
cx q[3],q[2];
rz(1.84207773208618) q[2];
sx q[2];
rz(5.18654361565644) q[2];
sx q[2];
rz(12.5659103155057) q[2];
rz(5.60962295532227) q[3];
sx q[3];
rz(5.3847703059488) q[3];
sx q[3];
rz(11.8571777105252) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.85072684288025) q[0];
sx q[0];
rz(4.45575872262055) q[0];
sx q[0];
rz(10.7283962726514) q[0];
rz(-0.297790080308914) q[1];
sx q[1];
rz(4.03715899785096) q[1];
sx q[1];
rz(12.272490477554) q[1];
cx q[1],q[0];
rz(1.03369081020355) q[0];
sx q[0];
rz(-2.14060386816924) q[0];
sx q[0];
rz(10.3305678129117) q[0];
rz(-2.00194644927979) q[2];
sx q[2];
rz(5.34629813035066) q[2];
sx q[2];
rz(5.23627135752841) q[2];
cx q[2],q[1];
rz(3.30505585670471) q[1];
sx q[1];
rz(2.46789643366868) q[1];
sx q[1];
rz(13.3188037633817) q[1];
rz(4.07892465591431) q[3];
sx q[3];
rz(6.66804805596406) q[3];
sx q[3];
rz(11.8402766942899) q[3];
cx q[3],q[2];
rz(-0.46875724196434) q[2];
sx q[2];
rz(4.57772043545777) q[2];
sx q[2];
rz(12.7871472597043) q[2];
rz(1.86864340305328) q[3];
sx q[3];
rz(4.88844552834565) q[3];
sx q[3];
rz(11.4181970119397) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.4918167591095) q[0];
sx q[0];
rz(0.621155412989207) q[0];
sx q[0];
rz(8.85028526782199) q[0];
rz(-0.542213976383209) q[1];
sx q[1];
rz(1.50339904625947) q[1];
sx q[1];
rz(6.67393205165073) q[1];
cx q[1],q[0];
rz(-0.465136557817459) q[0];
sx q[0];
rz(5.30445233185823) q[0];
sx q[0];
rz(9.59916793405219) q[0];
rz(1.62298774719238) q[2];
sx q[2];
rz(4.16858425934846) q[2];
sx q[2];
rz(11.7700724363248) q[2];
cx q[2],q[1];
rz(3.86061382293701) q[1];
sx q[1];
rz(3.55836323102052) q[1];
sx q[1];
rz(7.86762509345218) q[1];
rz(1.47467958927155) q[3];
sx q[3];
rz(6.65046253998811) q[3];
sx q[3];
rz(9.15099341272517) q[3];
cx q[3],q[2];
rz(1.6812402009964) q[2];
sx q[2];
rz(1.94430235226686) q[2];
sx q[2];
rz(11.9536199331205) q[2];
rz(2.15933060646057) q[3];
sx q[3];
rz(4.96862414677674) q[3];
sx q[3];
rz(13.0314740896146) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.149234443902969) q[0];
sx q[0];
rz(4.71639374096925) q[0];
sx q[0];
rz(9.36944311707422) q[0];
rz(1.55703592300415) q[1];
sx q[1];
rz(4.22959974606568) q[1];
sx q[1];
rz(6.72311351298496) q[1];
cx q[1],q[0];
rz(-1.20383679866791) q[0];
sx q[0];
rz(1.81696620781953) q[0];
sx q[0];
rz(7.38105342387363) q[0];
rz(0.20452044904232) q[2];
sx q[2];
rz(5.1946908553415) q[2];
sx q[2];
rz(8.74789986609622) q[2];
cx q[2],q[1];
rz(-1.4127368927002) q[1];
sx q[1];
rz(4.3970959504419) q[1];
sx q[1];
rz(13.5562920331876) q[1];
rz(1.85274708271027) q[3];
sx q[3];
rz(8.40226093133027) q[3];
sx q[3];
rz(7.86152539252445) q[3];
cx q[3],q[2];
rz(0.767919301986694) q[2];
sx q[2];
rz(4.93506375153596) q[2];
sx q[2];
rz(9.11276060938045) q[2];
rz(0.919655203819275) q[3];
sx q[3];
rz(4.91474524338777) q[3];
sx q[3];
rz(9.64095022379562) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.311740815639496) q[0];
sx q[0];
rz(2.15992018778855) q[0];
sx q[0];
rz(9.45290697588726) q[0];
rz(4.83713150024414) q[1];
sx q[1];
rz(1.96069339116151) q[1];
sx q[1];
rz(8.96528071760341) q[1];
cx q[1],q[0];
rz(1.05796790122986) q[0];
sx q[0];
rz(3.09593847219879) q[0];
sx q[0];
rz(10.0999386668126) q[0];
rz(-0.358768224716187) q[2];
sx q[2];
rz(3.9229166825586) q[2];
sx q[2];
rz(7.89748749732181) q[2];
cx q[2],q[1];
rz(2.50024366378784) q[1];
sx q[1];
rz(5.38796320756013) q[1];
sx q[1];
rz(7.66828129290744) q[1];
rz(4.57313442230225) q[3];
sx q[3];
rz(4.33365848858888) q[3];
sx q[3];
rz(8.26832065581485) q[3];
cx q[3],q[2];
rz(0.100998990237713) q[2];
sx q[2];
rz(4.93855610688264) q[2];
sx q[2];
rz(9.67994383572742) q[2];
rz(-0.986628711223602) q[3];
sx q[3];
rz(5.86846295197541) q[3];
sx q[3];
rz(10.8478759288709) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.76070177555084) q[0];
sx q[0];
rz(2.07168999512727) q[0];
sx q[0];
rz(8.08263728617831) q[0];
rz(-0.0019207150908187) q[1];
sx q[1];
rz(5.24058857758576) q[1];
sx q[1];
rz(5.23326728343173) q[1];
cx q[1],q[0];
rz(0.382978886365891) q[0];
sx q[0];
rz(2.88141238887841) q[0];
sx q[0];
rz(10.0364029169004) q[0];
rz(1.57200717926025) q[2];
sx q[2];
rz(4.8883537371927) q[2];
sx q[2];
rz(7.03544948100253) q[2];
cx q[2],q[1];
rz(8.39278030395508) q[1];
sx q[1];
rz(1.20326152642305) q[1];
sx q[1];
rz(6.04930756091281) q[1];
rz(-0.240137338638306) q[3];
sx q[3];
rz(3.88615080912644) q[3];
sx q[3];
rz(8.2575539112012) q[3];
cx q[3],q[2];
rz(0.929010152816772) q[2];
sx q[2];
rz(1.92761030991609) q[2];
sx q[2];
rz(13.0623671770017) q[2];
rz(1.46047329902649) q[3];
sx q[3];
rz(4.48335257371003) q[3];
sx q[3];
rz(8.23783383368655) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0542520508170128) q[0];
sx q[0];
rz(2.65340092976625) q[0];
sx q[0];
rz(8.14348444937869) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(0.651650846004486) q[1];
sx q[1];
rz(2.28545060952241) q[1];
sx q[1];
rz(9.83824727534457) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(3.10508394241333) q[2];
sx q[2];
rz(2.6513228734308) q[2];
sx q[2];
rz(5.79491469859287) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.913367092609406) q[3];
sx q[3];
rz(4.37521675427491) q[3];
sx q[3];
rz(16.3940233945768) q[3];
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
