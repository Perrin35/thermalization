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
rz(-0.705132007598877) q[0];
sx q[0];
rz(6.83505144913728) q[0];
sx q[0];
rz(9.44661363250717) q[0];
rz(-0.394373744726181) q[1];
sx q[1];
rz(4.60125759442384) q[1];
sx q[1];
rz(9.63969456254646) q[1];
cx q[1],q[0];
rz(-0.414799332618713) q[0];
sx q[0];
rz(2.7601294537359) q[0];
sx q[0];
rz(7.03857443331882) q[0];
rz(5.01073789596558) q[2];
sx q[2];
rz(3.62924447854096) q[2];
sx q[2];
rz(5.01144120692416) q[2];
cx q[2],q[1];
rz(-0.546751976013184) q[1];
sx q[1];
rz(1.12101236184175) q[1];
sx q[1];
rz(11.9064070939939) q[1];
rz(-0.488223701715469) q[3];
sx q[3];
rz(4.05752006371553) q[3];
sx q[3];
rz(10.5141130447309) q[3];
cx q[3],q[2];
rz(-0.731378376483917) q[2];
sx q[2];
rz(4.60089948971803) q[2];
sx q[2];
rz(13.1305756330411) q[2];
rz(-1.36518597602844) q[3];
sx q[3];
rz(3.59122148354585) q[3];
sx q[3];
rz(7.55243501662418) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.04417014122009) q[0];
sx q[0];
rz(1.92822697957093) q[0];
sx q[0];
rz(10.3527608871381) q[0];
rz(-1.16529750823975) q[1];
sx q[1];
rz(4.74492100079591) q[1];
sx q[1];
rz(11.6696197748105) q[1];
cx q[1],q[0];
rz(2.35515356063843) q[0];
sx q[0];
rz(5.54698243935639) q[0];
sx q[0];
rz(12.2720963716428) q[0];
rz(1.4019581079483) q[2];
sx q[2];
rz(1.64279201825196) q[2];
sx q[2];
rz(12.8333046197812) q[2];
cx q[2],q[1];
rz(1.08385443687439) q[1];
sx q[1];
rz(4.58855489094789) q[1];
sx q[1];
rz(11.7891094446103) q[1];
rz(-0.75150591135025) q[3];
sx q[3];
rz(1.25804010232026) q[3];
sx q[3];
rz(12.1137118101041) q[3];
cx q[3],q[2];
rz(-2.87598776817322) q[2];
sx q[2];
rz(5.74819818337495) q[2];
sx q[2];
rz(11.5262231588285) q[2];
rz(1.68637192249298) q[3];
sx q[3];
rz(1.29295948346192) q[3];
sx q[3];
rz(9.77953777312442) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.741376936435699) q[0];
sx q[0];
rz(2.56403604348237) q[0];
sx q[0];
rz(11.5381691217343) q[0];
rz(1.07854115962982) q[1];
sx q[1];
rz(2.57873472769792) q[1];
sx q[1];
rz(9.85991198419734) q[1];
cx q[1],q[0];
rz(3.80575656890869) q[0];
sx q[0];
rz(7.74648967583711) q[0];
sx q[0];
rz(10.6420739650647) q[0];
rz(1.82110941410065) q[2];
sx q[2];
rz(4.56969097455079) q[2];
sx q[2];
rz(6.66880891322299) q[2];
cx q[2],q[1];
rz(2.03554654121399) q[1];
sx q[1];
rz(4.54910078843171) q[1];
sx q[1];
rz(9.91307202576801) q[1];
rz(0.80743396282196) q[3];
sx q[3];
rz(4.9660713990503) q[3];
sx q[3];
rz(12.2201258897702) q[3];
cx q[3],q[2];
rz(0.533263981342316) q[2];
sx q[2];
rz(4.44190576870973) q[2];
sx q[2];
rz(6.58610055445834) q[2];
rz(1.32519245147705) q[3];
sx q[3];
rz(1.15851000149781) q[3];
sx q[3];
rz(9.33375286906167) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.94516968727112) q[0];
sx q[0];
rz(1.4099932034784) q[0];
sx q[0];
rz(10.3422280907552) q[0];
rz(-2.46871852874756) q[1];
sx q[1];
rz(5.1976901610666) q[1];
sx q[1];
rz(9.15990222095653) q[1];
cx q[1],q[0];
rz(-1.05658555030823) q[0];
sx q[0];
rz(4.68228998978669) q[0];
sx q[0];
rz(8.5132467508237) q[0];
rz(1.84899234771729) q[2];
sx q[2];
rz(1.28312662442262) q[2];
sx q[2];
rz(10.3904371023099) q[2];
cx q[2],q[1];
rz(2.53238224983215) q[1];
sx q[1];
rz(4.01144102414186) q[1];
sx q[1];
rz(12.0246274232785) q[1];
rz(-2.35142397880554) q[3];
sx q[3];
rz(4.12299427588517) q[3];
sx q[3];
rz(11.336348271362) q[3];
cx q[3],q[2];
rz(0.363105684518814) q[2];
sx q[2];
rz(3.6299097259813) q[2];
sx q[2];
rz(10.9898140191953) q[2];
rz(2.11450839042664) q[3];
sx q[3];
rz(3.88306459982926) q[3];
sx q[3];
rz(11.4649913072507) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.890013694763184) q[0];
sx q[0];
rz(3.00478342373902) q[0];
sx q[0];
rz(9.90350959300205) q[0];
rz(2.10846734046936) q[1];
sx q[1];
rz(8.45354286034638) q[1];
sx q[1];
rz(8.47212538718387) q[1];
cx q[1],q[0];
rz(0.272723078727722) q[0];
sx q[0];
rz(2.58167436917359) q[0];
sx q[0];
rz(8.07627389430209) q[0];
rz(-0.907180428504944) q[2];
sx q[2];
rz(0.510820778208323) q[2];
sx q[2];
rz(6.6520452260892) q[2];
cx q[2],q[1];
rz(0.885412335395813) q[1];
sx q[1];
rz(3.60586142738397) q[1];
sx q[1];
rz(6.83534548281833) q[1];
rz(1.65145456790924) q[3];
sx q[3];
rz(5.80379787285859) q[3];
sx q[3];
rz(11.9978813886563) q[3];
cx q[3],q[2];
rz(1.71973633766174) q[2];
sx q[2];
rz(6.64599791367585) q[2];
sx q[2];
rz(13.0969519376676) q[2];
rz(1.73557078838348) q[3];
sx q[3];
rz(5.15501085122163) q[3];
sx q[3];
rz(8.59310952424213) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.567531645298004) q[0];
sx q[0];
rz(1.6447789986902) q[0];
sx q[0];
rz(11.0497379064481) q[0];
rz(4.9780797958374) q[1];
sx q[1];
rz(4.49248734314973) q[1];
sx q[1];
rz(6.11061284541293) q[1];
cx q[1],q[0];
rz(4.14623022079468) q[0];
sx q[0];
rz(7.43755117257173) q[0];
sx q[0];
rz(11.0104897975843) q[0];
rz(-1.83038783073425) q[2];
sx q[2];
rz(4.26942959626252) q[2];
sx q[2];
rz(10.7834359168927) q[2];
cx q[2],q[1];
rz(3.1291618347168) q[1];
sx q[1];
rz(5.85504284699494) q[1];
sx q[1];
rz(9.68041405677005) q[1];
rz(-1.25949490070343) q[3];
sx q[3];
rz(0.654003055887767) q[3];
sx q[3];
rz(7.52103922366306) q[3];
cx q[3],q[2];
rz(-0.837957859039307) q[2];
sx q[2];
rz(4.59562447865541) q[2];
sx q[2];
rz(10.5514137506406) q[2];
rz(-0.782223284244537) q[3];
sx q[3];
rz(4.37707737286622) q[3];
sx q[3];
rz(10.7627677678983) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.288502663373947) q[0];
sx q[0];
rz(3.44809374411637) q[0];
sx q[0];
rz(10.0862537980001) q[0];
rz(5.32278966903687) q[1];
sx q[1];
rz(4.88211515744264) q[1];
sx q[1];
rz(13.3229698896329) q[1];
cx q[1],q[0];
rz(-6.62494230270386) q[0];
sx q[0];
rz(3.31986662943894) q[0];
sx q[0];
rz(6.83295390605136) q[0];
rz(-0.201621443033218) q[2];
sx q[2];
rz(5.11634353001649) q[2];
sx q[2];
rz(9.18555035292312) q[2];
cx q[2],q[1];
rz(3.50339889526367) q[1];
sx q[1];
rz(9.14531722863252) q[1];
sx q[1];
rz(11.9642584085385) q[1];
rz(0.0593341216444969) q[3];
sx q[3];
rz(1.37951281865174) q[3];
sx q[3];
rz(8.78292194604083) q[3];
cx q[3],q[2];
rz(-2.00442719459534) q[2];
sx q[2];
rz(6.09286990960176) q[2];
sx q[2];
rz(12.7608015298764) q[2];
rz(0.913131773471832) q[3];
sx q[3];
rz(4.89558580716188) q[3];
sx q[3];
rz(8.43936434983417) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.596543550491333) q[0];
sx q[0];
rz(2.52484348614747) q[0];
sx q[0];
rz(12.6330373048703) q[0];
rz(2.81702589988708) q[1];
sx q[1];
rz(4.77876702149446) q[1];
sx q[1];
rz(4.13047311305209) q[1];
cx q[1],q[0];
rz(2.75224041938782) q[0];
sx q[0];
rz(1.17124954064424) q[0];
sx q[0];
rz(10.6366208553235) q[0];
rz(-1.47739481925964) q[2];
sx q[2];
rz(4.7283161004358) q[2];
sx q[2];
rz(11.2956337690274) q[2];
cx q[2],q[1];
rz(2.73720717430115) q[1];
sx q[1];
rz(2.721524091559) q[1];
sx q[1];
rz(6.11893699168369) q[1];
rz(-0.819715797901154) q[3];
sx q[3];
rz(2.48206475575502) q[3];
sx q[3];
rz(10.2905382871549) q[3];
cx q[3],q[2];
rz(-1.46183347702026) q[2];
sx q[2];
rz(2.24488482077653) q[2];
sx q[2];
rz(12.9740040063779) q[2];
rz(-0.768610119819641) q[3];
sx q[3];
rz(4.9694369157129) q[3];
sx q[3];
rz(8.50087985991641) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.758936285972595) q[0];
sx q[0];
rz(4.98126080830628) q[0];
sx q[0];
rz(8.81557557582065) q[0];
rz(-3.23669743537903) q[1];
sx q[1];
rz(5.03118041356141) q[1];
sx q[1];
rz(10.2981551051061) q[1];
cx q[1],q[0];
rz(1.09268891811371) q[0];
sx q[0];
rz(2.94982189138467) q[0];
sx q[0];
rz(9.06345022319957) q[0];
rz(4.56486320495605) q[2];
sx q[2];
rz(4.57749822934205) q[2];
sx q[2];
rz(12.3019540071408) q[2];
cx q[2],q[1];
rz(-3.12102913856506) q[1];
sx q[1];
rz(3.82354071934754) q[1];
sx q[1];
rz(13.6374325513761) q[1];
rz(-2.01054763793945) q[3];
sx q[3];
rz(5.24506750901277) q[3];
sx q[3];
rz(9.03567314743205) q[3];
cx q[3],q[2];
rz(2.02188992500305) q[2];
sx q[2];
rz(7.07122054894502) q[2];
sx q[2];
rz(10.1070818662564) q[2];
rz(0.374263793230057) q[3];
sx q[3];
rz(1.57286098797853) q[3];
sx q[3];
rz(9.59168774484798) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.44362139701843) q[0];
sx q[0];
rz(4.5014022906595) q[0];
sx q[0];
rz(8.47181067465946) q[0];
rz(0.826418101787567) q[1];
sx q[1];
rz(3.88077143033082) q[1];
sx q[1];
rz(11.1698940753858) q[1];
cx q[1],q[0];
rz(0.291659265756607) q[0];
sx q[0];
rz(3.47591918905313) q[0];
sx q[0];
rz(10.7313800811689) q[0];
rz(1.94864046573639) q[2];
sx q[2];
rz(2.44777384598786) q[2];
sx q[2];
rz(9.99048224686786) q[2];
cx q[2],q[1];
rz(1.6203305721283) q[1];
sx q[1];
rz(4.42790296872193) q[1];
sx q[1];
rz(6.60744783877536) q[1];
rz(0.0568043813109398) q[3];
sx q[3];
rz(1.83265224297578) q[3];
sx q[3];
rz(8.23541567324802) q[3];
cx q[3],q[2];
rz(-2.67541265487671) q[2];
sx q[2];
rz(0.356234939890452) q[2];
sx q[2];
rz(9.58456811904117) q[2];
rz(-0.301884889602661) q[3];
sx q[3];
rz(5.35620632966096) q[3];
sx q[3];
rz(6.5704281091611) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-3.07821011543274) q[0];
sx q[0];
rz(2.2463255246454) q[0];
sx q[0];
rz(7.84219477175876) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-1.44443213939667) q[1];
sx q[1];
rz(1.44241383870179) q[1];
sx q[1];
rz(9.37084539084836) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.864669978618622) q[2];
sx q[2];
rz(5.16623106797273) q[2];
sx q[2];
rz(9.25395918487712) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.566833257675171) q[3];
sx q[3];
rz(4.36750045617158) q[3];
sx q[3];
rz(10.5778037071149) q[3];
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
