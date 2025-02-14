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
rz(0.604316890239716) q[0];
sx q[0];
rz(3.38316264946992) q[0];
sx q[0];
rz(9.09175050853893) q[0];
rz(-1.13554072380066) q[1];
sx q[1];
rz(3.96852132876451) q[1];
sx q[1];
rz(10.0687383770864) q[1];
cx q[1],q[0];
rz(0.0523833110928535) q[0];
sx q[0];
rz(4.10873827536637) q[0];
sx q[0];
rz(10.6643661022107) q[0];
rz(2.04333019256592) q[2];
sx q[2];
rz(2.80715495546395) q[2];
sx q[2];
rz(8.84485862254306) q[2];
cx q[2],q[1];
rz(-0.871306896209717) q[1];
sx q[1];
rz(4.47568491299684) q[1];
sx q[1];
rz(11.384137248985) q[1];
rz(0.412360101938248) q[3];
sx q[3];
rz(5.03318372567231) q[3];
sx q[3];
rz(9.29551990925475) q[3];
cx q[3],q[2];
rz(0.486462712287903) q[2];
sx q[2];
rz(3.60594475467736) q[2];
sx q[2];
rz(10.1655182003896) q[2];
rz(1.58985340595245) q[3];
sx q[3];
rz(3.85311028559739) q[3];
sx q[3];
rz(10.0175204634587) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.433166772127151) q[0];
sx q[0];
rz(4.27511504490907) q[0];
sx q[0];
rz(9.54174791871711) q[0];
rz(-0.504326939582825) q[1];
sx q[1];
rz(1.80289903481538) q[1];
sx q[1];
rz(9.70887622832462) q[1];
cx q[1],q[0];
rz(0.117825381457806) q[0];
sx q[0];
rz(4.62909379799897) q[0];
sx q[0];
rz(9.99154106377765) q[0];
rz(0.581677258014679) q[2];
sx q[2];
rz(5.42043391068513) q[2];
sx q[2];
rz(9.31134074031516) q[2];
cx q[2],q[1];
rz(0.717754662036896) q[1];
sx q[1];
rz(1.53621295292909) q[1];
sx q[1];
rz(9.54688376783534) q[1];
rz(0.733621895313263) q[3];
sx q[3];
rz(3.50808400114114) q[3];
sx q[3];
rz(8.73819390534564) q[3];
cx q[3],q[2];
rz(0.326101928949356) q[2];
sx q[2];
rz(4.02719870408113) q[2];
sx q[2];
rz(10.1855592489164) q[2];
rz(0.869593620300293) q[3];
sx q[3];
rz(5.01709369023378) q[3];
sx q[3];
rz(8.97168085574313) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.780754387378693) q[0];
sx q[0];
rz(4.26190081437165) q[0];
sx q[0];
rz(9.67694858311816) q[0];
rz(1.44226562976837) q[1];
sx q[1];
rz(2.18630543549592) q[1];
sx q[1];
rz(9.06851304172679) q[1];
cx q[1],q[0];
rz(0.399224549531937) q[0];
sx q[0];
rz(3.08168735553557) q[0];
sx q[0];
rz(9.32128630428716) q[0];
rz(-0.411150723695755) q[2];
sx q[2];
rz(4.20142653782899) q[2];
sx q[2];
rz(11.4069356679837) q[2];
cx q[2],q[1];
rz(0.608818471431732) q[1];
sx q[1];
rz(4.57007363637025) q[1];
sx q[1];
rz(10.0641909599225) q[1];
rz(0.00744354445487261) q[3];
sx q[3];
rz(3.39791637857492) q[3];
sx q[3];
rz(10.5074927568357) q[3];
cx q[3],q[2];
rz(1.11892712116241) q[2];
sx q[2];
rz(4.89291218121583) q[2];
sx q[2];
rz(11.0537911415021) q[2];
rz(0.00969836115837097) q[3];
sx q[3];
rz(5.05278745492036) q[3];
sx q[3];
rz(10.0461787938993) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.64647221565247) q[0];
sx q[0];
rz(3.68397394021089) q[0];
sx q[0];
rz(9.67553979753658) q[0];
rz(1.79509019851685) q[1];
sx q[1];
rz(5.75400462945039) q[1];
sx q[1];
rz(11.1231103897016) q[1];
cx q[1],q[0];
rz(-1.29372096061707) q[0];
sx q[0];
rz(4.92847982247407) q[0];
sx q[0];
rz(8.67306683062717) q[0];
rz(0.629486739635468) q[2];
sx q[2];
rz(2.51385900576646) q[2];
sx q[2];
rz(9.66772193311855) q[2];
cx q[2],q[1];
rz(1.55746746063232) q[1];
sx q[1];
rz(2.26557764609391) q[1];
sx q[1];
rz(7.13502285479709) q[1];
rz(-0.269142925739288) q[3];
sx q[3];
rz(1.82401040394837) q[3];
sx q[3];
rz(11.3879531383435) q[3];
cx q[3],q[2];
rz(-0.555633902549744) q[2];
sx q[2];
rz(5.17796793778474) q[2];
sx q[2];
rz(9.18146114646598) q[2];
rz(0.584687411785126) q[3];
sx q[3];
rz(3.71875378687913) q[3];
sx q[3];
rz(9.39664436354443) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.99341595172882) q[0];
sx q[0];
rz(3.50734093983705) q[0];
sx q[0];
rz(9.60433763860866) q[0];
rz(1.91636323928833) q[1];
sx q[1];
rz(4.5461175759607) q[1];
sx q[1];
rz(8.96652037500545) q[1];
cx q[1],q[0];
rz(0.970446348190308) q[0];
sx q[0];
rz(4.41685727437074) q[0];
sx q[0];
rz(9.34350211768552) q[0];
rz(0.266368359327316) q[2];
sx q[2];
rz(2.28122160037095) q[2];
sx q[2];
rz(11.3480119466703) q[2];
cx q[2],q[1];
rz(-1.04416644573212) q[1];
sx q[1];
rz(1.49617448647554) q[1];
sx q[1];
rz(9.55973712205097) q[1];
rz(0.700883746147156) q[3];
sx q[3];
rz(4.43220153649385) q[3];
sx q[3];
rz(10.2537137031476) q[3];
cx q[3],q[2];
rz(-1.75635838508606) q[2];
sx q[2];
rz(3.56591987808282) q[2];
sx q[2];
rz(10.3465620636861) q[2];
rz(1.890016913414) q[3];
sx q[3];
rz(4.87488904793794) q[3];
sx q[3];
rz(9.48644791393682) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0940304100513458) q[0];
sx q[0];
rz(4.05377629597718) q[0];
sx q[0];
rz(9.27611274122401) q[0];
rz(2.8246762752533) q[1];
sx q[1];
rz(5.80444851716096) q[1];
sx q[1];
rz(8.80794308184787) q[1];
cx q[1],q[0];
rz(1.51825225353241) q[0];
sx q[0];
rz(3.04193110962445) q[0];
sx q[0];
rz(8.7600250005643) q[0];
rz(-0.440321147441864) q[2];
sx q[2];
rz(4.38402751286561) q[2];
sx q[2];
rz(11.6963693857114) q[2];
cx q[2],q[1];
rz(1.17795360088348) q[1];
sx q[1];
rz(5.65278163750703) q[1];
sx q[1];
rz(9.85823977588817) q[1];
rz(-0.455801397562027) q[3];
sx q[3];
rz(3.71304968197877) q[3];
sx q[3];
rz(9.14018697141811) q[3];
cx q[3],q[2];
rz(0.864691197872162) q[2];
sx q[2];
rz(4.57029763062532) q[2];
sx q[2];
rz(9.43310423231825) q[2];
rz(0.626380383968353) q[3];
sx q[3];
rz(3.85758659442002) q[3];
sx q[3];
rz(9.42532869860671) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.3061226606369) q[0];
sx q[0];
rz(5.13130393822724) q[0];
sx q[0];
rz(9.90685949324771) q[0];
rz(-0.0291906334459782) q[1];
sx q[1];
rz(1.84552851517732) q[1];
sx q[1];
rz(10.1888203978459) q[1];
cx q[1],q[0];
rz(0.692248344421387) q[0];
sx q[0];
rz(1.63817611535127) q[0];
sx q[0];
rz(9.79601824878856) q[0];
rz(0.0326235517859459) q[2];
sx q[2];
rz(2.4815686066919) q[2];
sx q[2];
rz(11.7775203943174) q[2];
cx q[2],q[1];
rz(-0.936385691165924) q[1];
sx q[1];
rz(4.11870679457719) q[1];
sx q[1];
rz(11.200632429115) q[1];
rz(0.0188347119837999) q[3];
sx q[3];
rz(5.88220134575898) q[3];
sx q[3];
rz(9.74141854643031) q[3];
cx q[3],q[2];
rz(0.463341057300568) q[2];
sx q[2];
rz(4.96387139161164) q[2];
sx q[2];
rz(9.90894440411731) q[2];
rz(-0.682289958000183) q[3];
sx q[3];
rz(3.79569348891313) q[3];
sx q[3];
rz(9.29003762303993) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0844556465744972) q[0];
sx q[0];
rz(3.91939309437806) q[0];
sx q[0];
rz(8.13301119803592) q[0];
rz(-0.465033978223801) q[1];
sx q[1];
rz(3.66042307217652) q[1];
sx q[1];
rz(9.73194188474818) q[1];
cx q[1],q[0];
rz(0.391246736049652) q[0];
sx q[0];
rz(5.21747437317903) q[0];
sx q[0];
rz(11.4822077512662) q[0];
rz(0.620218992233276) q[2];
sx q[2];
rz(2.38313159544999) q[2];
sx q[2];
rz(11.1454479455869) q[2];
cx q[2],q[1];
rz(-1.39138448238373) q[1];
sx q[1];
rz(2.54738745291764) q[1];
sx q[1];
rz(9.7069341301839) q[1];
rz(0.0203367881476879) q[3];
sx q[3];
rz(5.33738747437532) q[3];
sx q[3];
rz(9.96381059884235) q[3];
cx q[3],q[2];
rz(0.180535987019539) q[2];
sx q[2];
rz(2.70132160385186) q[2];
sx q[2];
rz(8.70467886923953) q[2];
rz(0.850472033023834) q[3];
sx q[3];
rz(4.30840733845765) q[3];
sx q[3];
rz(7.69478175639316) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.20486031472683) q[0];
sx q[0];
rz(3.31008038123185) q[0];
sx q[0];
rz(10.1329244732778) q[0];
rz(1.5180960893631) q[1];
sx q[1];
rz(4.0797462781244) q[1];
sx q[1];
rz(9.06858233212634) q[1];
cx q[1],q[0];
rz(-0.307638764381409) q[0];
sx q[0];
rz(3.81265846093232) q[0];
sx q[0];
rz(10.2789397597234) q[0];
rz(-0.0296102799475193) q[2];
sx q[2];
rz(5.3300680239969) q[2];
sx q[2];
rz(9.06707180141612) q[2];
cx q[2],q[1];
rz(0.569363832473755) q[1];
sx q[1];
rz(3.79112103779847) q[1];
sx q[1];
rz(7.80429408549472) q[1];
rz(-0.709077656269073) q[3];
sx q[3];
rz(3.72522279818589) q[3];
sx q[3];
rz(10.2511310338895) q[3];
cx q[3],q[2];
rz(1.09229624271393) q[2];
sx q[2];
rz(5.47587433655793) q[2];
sx q[2];
rz(10.3105699181478) q[2];
rz(2.20573353767395) q[3];
sx q[3];
rz(4.32212344010408) q[3];
sx q[3];
rz(9.6460535287778) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.437792718410492) q[0];
sx q[0];
rz(3.09428230871493) q[0];
sx q[0];
rz(10.8742930650632) q[0];
rz(-1.02251470088959) q[1];
sx q[1];
rz(2.65785637696321) q[1];
sx q[1];
rz(11.1170695781629) q[1];
cx q[1],q[0];
rz(-0.677748620510101) q[0];
sx q[0];
rz(2.72632563312585) q[0];
sx q[0];
rz(9.30144273339912) q[0];
rz(0.804331243038177) q[2];
sx q[2];
rz(3.70547828276689) q[2];
sx q[2];
rz(10.9248908519666) q[2];
cx q[2],q[1];
rz(0.305617988109589) q[1];
sx q[1];
rz(3.73753634293611) q[1];
sx q[1];
rz(9.8093411385934) q[1];
rz(-0.265492081642151) q[3];
sx q[3];
rz(1.42783001263673) q[3];
sx q[3];
rz(13.1056809186856) q[3];
cx q[3],q[2];
rz(-0.255410432815552) q[2];
sx q[2];
rz(4.33151260216767) q[2];
sx q[2];
rz(10.0996463656346) q[2];
rz(2.17009615898132) q[3];
sx q[3];
rz(2.43923047383363) q[3];
sx q[3];
rz(7.93318996428653) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.199344247579575) q[0];
sx q[0];
rz(3.85540422995622) q[0];
sx q[0];
rz(9.46310890316173) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(0.157507732510567) q[1];
sx q[1];
rz(3.73523786862428) q[1];
sx q[1];
rz(7.52566657065555) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-2*pi/9) q[2];
sx q[2];
rz(4.87826720078523) q[2];
sx q[2];
rz(10.3820282578389) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-1.07157647609711) q[3];
sx q[3];
rz(4.2518951018625) q[3];
sx q[3];
rz(11.5847074747007) q[3];
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
