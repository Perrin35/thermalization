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
rz(1.00320494174957) q[0];
sx q[0];
rz(1.94110742409761) q[0];
sx q[0];
rz(8.39841256140872) q[0];
rz(0.298445552587509) q[1];
sx q[1];
rz(4.65000000794465) q[1];
sx q[1];
rz(10.4273726701657) q[1];
cx q[1],q[0];
rz(-0.745131969451904) q[0];
sx q[0];
rz(3.34667439957196) q[0];
sx q[0];
rz(9.53701422958776) q[0];
rz(-0.907881259918213) q[2];
sx q[2];
rz(1.61432448227937) q[2];
sx q[2];
rz(11.2211041212003) q[2];
cx q[2],q[1];
rz(0.943843126296997) q[1];
sx q[1];
rz(1.70438662369783) q[1];
sx q[1];
rz(11.3859109640042) q[1];
rz(-2.53462362289429) q[3];
sx q[3];
rz(4.16983249981935) q[3];
sx q[3];
rz(11.7822496652524) q[3];
cx q[3],q[2];
rz(0.588743209838867) q[2];
sx q[2];
rz(2.14079490502412) q[2];
sx q[2];
rz(12.6014783143918) q[2];
rz(1.2037752866745) q[3];
sx q[3];
rz(4.96098450024659) q[3];
sx q[3];
rz(9.70903647541209) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.069000706076622) q[0];
sx q[0];
rz(3.5836537202173) q[0];
sx q[0];
rz(10.4289850950162) q[0];
rz(-0.122411571443081) q[1];
sx q[1];
rz(4.60917690594728) q[1];
sx q[1];
rz(8.76769814490482) q[1];
cx q[1],q[0];
rz(1.94170379638672) q[0];
sx q[0];
rz(3.53797841270501) q[0];
sx q[0];
rz(8.06740329264804) q[0];
rz(-1.73209249973297) q[2];
sx q[2];
rz(4.9862928708368) q[2];
sx q[2];
rz(11.5581705331723) q[2];
cx q[2],q[1];
rz(-2.70656943321228) q[1];
sx q[1];
rz(3.98865977128083) q[1];
sx q[1];
rz(13.2212524175565) q[1];
rz(0.244796276092529) q[3];
sx q[3];
rz(3.94362202485139) q[3];
sx q[3];
rz(10.5001240730207) q[3];
cx q[3],q[2];
rz(1.58480262756348) q[2];
sx q[2];
rz(5.70527234871919) q[2];
sx q[2];
rz(11.9744879960935) q[2];
rz(1.66517341136932) q[3];
sx q[3];
rz(4.74905279477174) q[3];
sx q[3];
rz(11.6927230119626) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.00374352931976) q[0];
sx q[0];
rz(3.35054487188394) q[0];
sx q[0];
rz(10.6393999814908) q[0];
rz(2.18517231941223) q[1];
sx q[1];
rz(1.95018699963624) q[1];
sx q[1];
rz(8.2531536579053) q[1];
cx q[1],q[0];
rz(1.71922433376312) q[0];
sx q[0];
rz(3.77243599494035) q[0];
sx q[0];
rz(10.1331901907842) q[0];
rz(0.700525224208832) q[2];
sx q[2];
rz(4.8702051957422) q[2];
sx q[2];
rz(11.1036567449491) q[2];
cx q[2],q[1];
rz(0.708327889442444) q[1];
sx q[1];
rz(4.37557211716706) q[1];
sx q[1];
rz(10.3958621978681) q[1];
rz(0.0449186712503433) q[3];
sx q[3];
rz(4.62262574036653) q[3];
sx q[3];
rz(10.3843855023305) q[3];
cx q[3],q[2];
rz(0.200229153037071) q[2];
sx q[2];
rz(4.1919221003824) q[2];
sx q[2];
rz(10.8024558782498) q[2];
rz(0.406698346138) q[3];
sx q[3];
rz(3.79087957938249) q[3];
sx q[3];
rz(12.10586140155) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.35206258296967) q[0];
sx q[0];
rz(5.02383402188356) q[0];
sx q[0];
rz(11.3258964776914) q[0];
rz(-0.718183040618896) q[1];
sx q[1];
rz(4.9572390635782) q[1];
sx q[1];
rz(9.18582715689346) q[1];
cx q[1],q[0];
rz(1.09319400787354) q[0];
sx q[0];
rz(2.82805475791032) q[0];
sx q[0];
rz(8.25623450278446) q[0];
rz(-0.85748302936554) q[2];
sx q[2];
rz(5.33326927025849) q[2];
sx q[2];
rz(11.9560184240262) q[2];
cx q[2],q[1];
rz(-0.0935557782649994) q[1];
sx q[1];
rz(2.84479201038415) q[1];
sx q[1];
rz(12.059304213516) q[1];
rz(1.03624176979065) q[3];
sx q[3];
rz(4.82814195950563) q[3];
sx q[3];
rz(8.75621221064731) q[3];
cx q[3],q[2];
rz(-0.516039848327637) q[2];
sx q[2];
rz(4.8896713574701) q[2];
sx q[2];
rz(7.62612590789005) q[2];
rz(-0.881377339363098) q[3];
sx q[3];
rz(3.31048241456086) q[3];
sx q[3];
rz(10.2057708859365) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.41370224952698) q[0];
sx q[0];
rz(3.38562202652032) q[0];
sx q[0];
rz(10.9090820312421) q[0];
rz(-0.653857469558716) q[1];
sx q[1];
rz(1.52127185662324) q[1];
sx q[1];
rz(9.56271559595271) q[1];
cx q[1],q[0];
rz(1.11780953407288) q[0];
sx q[0];
rz(4.35950961907441) q[0];
sx q[0];
rz(10.8369868755262) q[0];
rz(-1.89100778102875) q[2];
sx q[2];
rz(3.81971278985078) q[2];
sx q[2];
rz(7.72434923648044) q[2];
cx q[2],q[1];
rz(-1.55224621295929) q[1];
sx q[1];
rz(4.81033566792543) q[1];
sx q[1];
rz(8.99330890773937) q[1];
rz(0.258853673934937) q[3];
sx q[3];
rz(4.76646104653413) q[3];
sx q[3];
rz(8.11720738410159) q[3];
cx q[3],q[2];
rz(0.589789271354675) q[2];
sx q[2];
rz(1.17295959790284) q[2];
sx q[2];
rz(9.51196239738866) q[2];
rz(-0.114635564386845) q[3];
sx q[3];
rz(3.95670291979844) q[3];
sx q[3];
rz(9.00747806429073) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.365994155406952) q[0];
sx q[0];
rz(4.74370745022828) q[0];
sx q[0];
rz(9.90974596738025) q[0];
rz(-0.943517744541168) q[1];
sx q[1];
rz(5.47019210656221) q[1];
sx q[1];
rz(10.643898344032) q[1];
cx q[1],q[0];
rz(0.232692256569862) q[0];
sx q[0];
rz(4.61498871644075) q[0];
sx q[0];
rz(10.1470138192098) q[0];
rz(1.58713829517365) q[2];
sx q[2];
rz(4.39040103753144) q[2];
sx q[2];
rz(13.6886257886808) q[2];
cx q[2],q[1];
rz(2.96129417419434) q[1];
sx q[1];
rz(4.83935943444306) q[1];
sx q[1];
rz(10.9190039396207) q[1];
rz(0.504849016666412) q[3];
sx q[3];
rz(5.40598407586152) q[3];
sx q[3];
rz(10.6170975923459) q[3];
cx q[3],q[2];
rz(-0.901482701301575) q[2];
sx q[2];
rz(4.20204010804231) q[2];
sx q[2];
rz(10.3660754322927) q[2];
rz(0.246020406484604) q[3];
sx q[3];
rz(1.49642077286775) q[3];
sx q[3];
rz(11.2062098741452) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0980023816227913) q[0];
sx q[0];
rz(5.18134728272492) q[0];
sx q[0];
rz(11.4023285865705) q[0];
rz(1.35804915428162) q[1];
sx q[1];
rz(4.0081349333101) q[1];
sx q[1];
rz(8.03174791335269) q[1];
cx q[1],q[0];
rz(1.00608015060425) q[0];
sx q[0];
rz(3.52515945036943) q[0];
sx q[0];
rz(9.6236713886182) q[0];
rz(-1.12474942207336) q[2];
sx q[2];
rz(2.54526332219178) q[2];
sx q[2];
rz(10.568732237808) q[2];
cx q[2],q[1];
rz(-1.60266149044037) q[1];
sx q[1];
rz(1.67625132401521) q[1];
sx q[1];
rz(11.4672996759336) q[1];
rz(1.79482531547546) q[3];
sx q[3];
rz(4.04600170453126) q[3];
sx q[3];
rz(10.6457102060239) q[3];
cx q[3],q[2];
rz(0.670365750789642) q[2];
sx q[2];
rz(4.85500017006929) q[2];
sx q[2];
rz(9.12582421898052) q[2];
rz(0.405458837747574) q[3];
sx q[3];
rz(2.77512276371057) q[3];
sx q[3];
rz(8.86045954226657) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.42951714992523) q[0];
sx q[0];
rz(2.86866644223268) q[0];
sx q[0];
rz(11.7770008802335) q[0];
rz(0.404965966939926) q[1];
sx q[1];
rz(4.49234214623506) q[1];
sx q[1];
rz(9.91973829864665) q[1];
cx q[1],q[0];
rz(-0.698236584663391) q[0];
sx q[0];
rz(3.91614362795884) q[0];
sx q[0];
rz(11.859663462631) q[0];
rz(1.27068603038788) q[2];
sx q[2];
rz(4.16296652157838) q[2];
sx q[2];
rz(13.5453042745511) q[2];
cx q[2],q[1];
rz(-0.725157618522644) q[1];
sx q[1];
rz(3.95924565394456) q[1];
sx q[1];
rz(10.6595767497937) q[1];
rz(0.734330952167511) q[3];
sx q[3];
rz(4.30691006978089) q[3];
sx q[3];
rz(8.74490789174243) q[3];
cx q[3],q[2];
rz(2.16275453567505) q[2];
sx q[2];
rz(1.91637221177156) q[2];
sx q[2];
rz(9.16952056287929) q[2];
rz(-0.0349595472216606) q[3];
sx q[3];
rz(4.66495898564393) q[3];
sx q[3];
rz(9.1785507261674) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.736647665500641) q[0];
sx q[0];
rz(5.22037163575227) q[0];
sx q[0];
rz(10.1348086357038) q[0];
rz(1.00114774703979) q[1];
sx q[1];
rz(4.03877958853776) q[1];
sx q[1];
rz(8.80387691258594) q[1];
cx q[1],q[0];
rz(1.37085211277008) q[0];
sx q[0];
rz(1.93198791344697) q[0];
sx q[0];
rz(11.0600906371991) q[0];
rz(-2.22775197029114) q[2];
sx q[2];
rz(1.95994201500947) q[2];
sx q[2];
rz(10.2126592755239) q[2];
cx q[2],q[1];
rz(0.945407509803772) q[1];
sx q[1];
rz(5.05923751195008) q[1];
sx q[1];
rz(9.31065001933976) q[1];
rz(-0.277554661035538) q[3];
sx q[3];
rz(3.955091329413) q[3];
sx q[3];
rz(9.39438218473598) q[3];
cx q[3],q[2];
rz(1.33224773406982) q[2];
sx q[2];
rz(1.52319088776643) q[2];
sx q[2];
rz(12.2771811246793) q[2];
rz(1.11221635341644) q[3];
sx q[3];
rz(3.43664613564546) q[3];
sx q[3];
rz(10.4451642990033) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.66388285160065) q[0];
sx q[0];
rz(1.79986706574494) q[0];
sx q[0];
rz(10.3515103816907) q[0];
rz(-1.1070739030838) q[1];
sx q[1];
rz(4.65534714062745) q[1];
sx q[1];
rz(9.99277094601795) q[1];
cx q[1],q[0];
rz(3.09997892379761) q[0];
sx q[0];
rz(5.02844408352906) q[0];
sx q[0];
rz(8.48810610770389) q[0];
rz(0.263754367828369) q[2];
sx q[2];
rz(4.21538785298402) q[2];
sx q[2];
rz(8.90121928452655) q[2];
cx q[2],q[1];
rz(0.895333290100098) q[1];
sx q[1];
rz(3.8833745439821) q[1];
sx q[1];
rz(9.98088405131503) q[1];
rz(0.915402114391327) q[3];
sx q[3];
rz(3.89565369685227) q[3];
sx q[3];
rz(9.47677381559416) q[3];
cx q[3],q[2];
rz(0.0522814355790615) q[2];
sx q[2];
rz(3.95189711649949) q[2];
sx q[2];
rz(7.98616180419132) q[2];
rz(-0.296640515327454) q[3];
sx q[3];
rz(3.48319280345971) q[3];
sx q[3];
rz(8.97773379682704) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.603070020675659) q[0];
sx q[0];
rz(2.07146159012849) q[0];
sx q[0];
rz(9.97107819317981) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(1.36645436286926) q[1];
sx q[1];
rz(6.20397058327729) q[1];
sx q[1];
rz(8.88035503625079) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(1.04713845252991) q[2];
sx q[2];
rz(3.81895241339738) q[2];
sx q[2];
rz(9.82401714324161) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.868114054203033) q[3];
sx q[3];
rz(5.01479974587495) q[3];
sx q[3];
rz(11.095348930351) q[3];
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
