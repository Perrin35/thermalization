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
rz(1.16601002216339) q[0];
sx q[0];
rz(2.26288852294023) q[0];
sx q[0];
rz(9.22876968084976) q[0];
rz(1.64855432510376) q[1];
sx q[1];
rz(4.07984075148637) q[1];
sx q[1];
rz(8.5473015665929) q[1];
cx q[1],q[0];
rz(3.39000344276428) q[0];
sx q[0];
rz(2.02010038693482) q[0];
sx q[0];
rz(11.2346510648648) q[0];
rz(-1.14669418334961) q[2];
sx q[2];
rz(4.39229586918885) q[2];
sx q[2];
rz(5.44170019625827) q[2];
cx q[2],q[1];
rz(-0.514982104301453) q[1];
sx q[1];
rz(4.26879993279512) q[1];
sx q[1];
rz(10.6969796180646) q[1];
rz(-0.438745021820068) q[3];
sx q[3];
rz(2.62257036765153) q[3];
sx q[3];
rz(10.6528796911161) q[3];
cx q[3],q[2];
rz(-0.332855463027954) q[2];
sx q[2];
rz(3.84940460522706) q[2];
sx q[2];
rz(12.5291781186978) q[2];
rz(0.913413643836975) q[3];
sx q[3];
rz(5.30089846451814) q[3];
sx q[3];
rz(7.79222545622989) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.344412416219711) q[0];
sx q[0];
rz(5.20887413819367) q[0];
sx q[0];
rz(9.20260786115333) q[0];
rz(3.33527946472168) q[1];
sx q[1];
rz(1.2682468016916) q[1];
sx q[1];
rz(9.16491675972148) q[1];
cx q[1],q[0];
rz(-1.3701137304306) q[0];
sx q[0];
rz(3.65538361867005) q[0];
sx q[0];
rz(9.93618944882556) q[0];
rz(-0.17879793047905) q[2];
sx q[2];
rz(4.8185135443979) q[2];
sx q[2];
rz(8.97896952032253) q[2];
cx q[2],q[1];
rz(10.1060237884521) q[1];
sx q[1];
rz(0.282398613291331) q[1];
sx q[1];
rz(10.2780950426976) q[1];
rz(0.81158322095871) q[3];
sx q[3];
rz(4.34622100194032) q[3];
sx q[3];
rz(9.42741229281902) q[3];
cx q[3],q[2];
rz(-0.766931772232056) q[2];
sx q[2];
rz(5.70034018357331) q[2];
sx q[2];
rz(9.58154957591697) q[2];
rz(0.805927038192749) q[3];
sx q[3];
rz(5.19278040726716) q[3];
sx q[3];
rz(9.94960907696887) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.904615163803101) q[0];
sx q[0];
rz(3.02930625726516) q[0];
sx q[0];
rz(11.3802325487058) q[0];
rz(-3.0778431892395) q[1];
sx q[1];
rz(4.66760888894136) q[1];
sx q[1];
rz(12.9788522481839) q[1];
cx q[1],q[0];
rz(0.291523993015289) q[0];
sx q[0];
rz(4.2253818829828) q[0];
sx q[0];
rz(10.0202048778455) q[0];
rz(1.72778761386871) q[2];
sx q[2];
rz(4.03848204215104) q[2];
sx q[2];
rz(11.0651110172193) q[2];
cx q[2],q[1];
rz(-0.361299157142639) q[1];
sx q[1];
rz(0.344432028132029) q[1];
sx q[1];
rz(13.4983920812528) q[1];
rz(2.15229773521423) q[3];
sx q[3];
rz(5.2921680529886) q[3];
sx q[3];
rz(10.6938934087674) q[3];
cx q[3],q[2];
rz(0.568466424942017) q[2];
sx q[2];
rz(3.58193740447099) q[2];
sx q[2];
rz(8.25863680838748) q[2];
rz(2.3447573184967) q[3];
sx q[3];
rz(4.91397956212098) q[3];
sx q[3];
rz(13.3476457357328) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.91549384593964) q[0];
sx q[0];
rz(4.77464929421479) q[0];
sx q[0];
rz(8.01139733790561) q[0];
rz(0.498902916908264) q[1];
sx q[1];
rz(4.52467289765412) q[1];
sx q[1];
rz(7.57705721854373) q[1];
cx q[1],q[0];
rz(-2.50996112823486) q[0];
sx q[0];
rz(-0.0892797390646258) q[0];
sx q[0];
rz(6.58882043360873) q[0];
rz(0.220458954572678) q[2];
sx q[2];
rz(3.47659963567788) q[2];
sx q[2];
rz(9.2230096667926) q[2];
cx q[2],q[1];
rz(-0.419997572898865) q[1];
sx q[1];
rz(4.65037957032258) q[1];
sx q[1];
rz(9.42796000804893) q[1];
rz(-2.31485152244568) q[3];
sx q[3];
rz(6.37477246125276) q[3];
sx q[3];
rz(10.5303562641065) q[3];
cx q[3],q[2];
rz(-0.529001951217651) q[2];
sx q[2];
rz(5.60219994385774) q[2];
sx q[2];
rz(12.7361814737241) q[2];
rz(0.428101927042007) q[3];
sx q[3];
rz(4.57110348542268) q[3];
sx q[3];
rz(9.58515969514056) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.47996461391449) q[0];
sx q[0];
rz(3.5339014848047) q[0];
sx q[0];
rz(11.7309913396756) q[0];
rz(-0.639404892921448) q[1];
sx q[1];
rz(3.80659315188462) q[1];
sx q[1];
rz(8.33379051684543) q[1];
cx q[1],q[0];
rz(3.3661892414093) q[0];
sx q[0];
rz(2.44516548712785) q[0];
sx q[0];
rz(7.87321350573703) q[0];
rz(-0.1674974411726) q[2];
sx q[2];
rz(7.55821529229219) q[2];
sx q[2];
rz(5.70492479800388) q[2];
cx q[2],q[1];
rz(-0.180269464850426) q[1];
sx q[1];
rz(1.1774447282129) q[1];
sx q[1];
rz(8.74534723757907) q[1];
rz(-1.897301197052) q[3];
sx q[3];
rz(4.21241894562776) q[3];
sx q[3];
rz(7.24477884768649) q[3];
cx q[3],q[2];
rz(-1.36844801902771) q[2];
sx q[2];
rz(1.1870563348108) q[2];
sx q[2];
rz(6.27914950846835) q[2];
rz(-1.07673799991608) q[3];
sx q[3];
rz(3.83919063408906) q[3];
sx q[3];
rz(12.7532746553342) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.31835973262787) q[0];
sx q[0];
rz(4.9218671639734) q[0];
sx q[0];
rz(8.79113451241657) q[0];
rz(3.54273581504822) q[1];
sx q[1];
rz(2.43234303792054) q[1];
sx q[1];
rz(6.97540972232028) q[1];
cx q[1],q[0];
rz(-0.617182314395905) q[0];
sx q[0];
rz(8.26470437844331) q[0];
sx q[0];
rz(10.1941528081815) q[0];
rz(-0.801815569400787) q[2];
sx q[2];
rz(4.48723903496797) q[2];
sx q[2];
rz(10.4482414483945) q[2];
cx q[2],q[1];
rz(-0.118973806500435) q[1];
sx q[1];
rz(4.46336451371247) q[1];
sx q[1];
rz(8.83812824486896) q[1];
rz(3.29254293441772) q[3];
sx q[3];
rz(4.52900162537629) q[3];
sx q[3];
rz(11.6517565011899) q[3];
cx q[3],q[2];
rz(2.93585181236267) q[2];
sx q[2];
rz(2.35710993607576) q[2];
sx q[2];
rz(8.18395493029758) q[2];
rz(1.95670175552368) q[3];
sx q[3];
rz(7.58457198937471) q[3];
sx q[3];
rz(14.1631822347562) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.0095089673996) q[0];
sx q[0];
rz(4.44460216363008) q[0];
sx q[0];
rz(10.7955552101056) q[0];
rz(-0.41807433962822) q[1];
sx q[1];
rz(4.62416795094544) q[1];
sx q[1];
rz(12.0797006845395) q[1];
cx q[1],q[0];
rz(0.80856990814209) q[0];
sx q[0];
rz(1.47367027600343) q[0];
sx q[0];
rz(8.86445114611789) q[0];
rz(-0.569966912269592) q[2];
sx q[2];
rz(2.51161328156526) q[2];
sx q[2];
rz(11.7153024434964) q[2];
cx q[2],q[1];
rz(-4.45188999176025) q[1];
sx q[1];
rz(6.01321354706819) q[1];
sx q[1];
rz(7.00535628794833) q[1];
rz(1.66310334205627) q[3];
sx q[3];
rz(0.47180167038972) q[3];
sx q[3];
rz(10.2863252520482) q[3];
cx q[3],q[2];
rz(0.7024986743927) q[2];
sx q[2];
rz(7.33217540581758) q[2];
sx q[2];
rz(12.7635669469754) q[2];
rz(0.509447336196899) q[3];
sx q[3];
rz(3.40224158962304) q[3];
sx q[3];
rz(8.59640929698154) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.86151230335236) q[0];
sx q[0];
rz(1.07416382630403) q[0];
sx q[0];
rz(8.0582941532056) q[0];
rz(-0.816614329814911) q[1];
sx q[1];
rz(2.05207029183442) q[1];
sx q[1];
rz(9.1420820414941) q[1];
cx q[1],q[0];
rz(0.824741005897522) q[0];
sx q[0];
rz(6.90877905686433) q[0];
sx q[0];
rz(9.73549581169292) q[0];
rz(1.85778915882111) q[2];
sx q[2];
rz(2.81761348445947) q[2];
sx q[2];
rz(8.50469300746127) q[2];
cx q[2],q[1];
rz(1.37794387340546) q[1];
sx q[1];
rz(6.99211970170076) q[1];
sx q[1];
rz(9.71305582522556) q[1];
rz(0.950403273105621) q[3];
sx q[3];
rz(3.75638953049714) q[3];
sx q[3];
rz(6.71841928958102) q[3];
cx q[3],q[2];
rz(0.0351598113775253) q[2];
sx q[2];
rz(4.14754191239411) q[2];
sx q[2];
rz(10.4079138398091) q[2];
rz(1.27849173545837) q[3];
sx q[3];
rz(5.35800662835176) q[3];
sx q[3];
rz(11.6374830961148) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.660166919231415) q[0];
sx q[0];
rz(4.93353405793244) q[0];
sx q[0];
rz(9.73804477452441) q[0];
rz(-3.66632127761841) q[1];
sx q[1];
rz(2.87867614825303) q[1];
sx q[1];
rz(11.3140365838925) q[1];
cx q[1],q[0];
rz(1.45626294612885) q[0];
sx q[0];
rz(2.09768167336518) q[0];
sx q[0];
rz(12.4267868757169) q[0];
rz(-3.94898772239685) q[2];
sx q[2];
rz(5.15546027024324) q[2];
sx q[2];
rz(9.32571141271993) q[2];
cx q[2],q[1];
rz(-2.55726671218872) q[1];
sx q[1];
rz(5.24460306962068) q[1];
sx q[1];
rz(13.5082454442899) q[1];
rz(1.65212035179138) q[3];
sx q[3];
rz(2.3235851248079) q[3];
sx q[3];
rz(9.9449649810712) q[3];
cx q[3],q[2];
rz(0.199204429984093) q[2];
sx q[2];
rz(5.65825262864167) q[2];
sx q[2];
rz(13.9478335142057) q[2];
rz(2.84046101570129) q[3];
sx q[3];
rz(4.12424442370469) q[3];
sx q[3];
rz(9.04376841186687) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.234066009521484) q[0];
sx q[0];
rz(2.13876870472962) q[0];
sx q[0];
rz(10.7722856759946) q[0];
rz(0.884959101676941) q[1];
sx q[1];
rz(5.6575335582071) q[1];
sx q[1];
rz(7.68166825770541) q[1];
cx q[1],q[0];
rz(0.325377255678177) q[0];
sx q[0];
rz(2.76698157389695) q[0];
sx q[0];
rz(9.39561758040591) q[0];
rz(0.705300688743591) q[2];
sx q[2];
rz(4.5719532092386) q[2];
sx q[2];
rz(12.8375499010007) q[2];
cx q[2],q[1];
rz(-0.815128862857819) q[1];
sx q[1];
rz(3.23020336230332) q[1];
sx q[1];
rz(13.7024445295255) q[1];
rz(0.353843778371811) q[3];
sx q[3];
rz(3.73305174906785) q[3];
sx q[3];
rz(6.50537226199313) q[3];
cx q[3],q[2];
rz(-3.19478297233582) q[2];
sx q[2];
rz(1.32064274151857) q[2];
sx q[2];
rz(10.3706260085027) q[2];
rz(-0.690190672874451) q[3];
sx q[3];
rz(4.3642282803827) q[3];
sx q[3];
rz(7.58420488833591) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.696619033813477) q[0];
sx q[0];
rz(7.14627710183198) q[0];
sx q[0];
rz(9.38193553908869) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(0.0892455950379372) q[1];
sx q[1];
rz(4.10239651997621) q[1];
sx q[1];
rz(10.9088883161466) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-2.71445274353027) q[2];
sx q[2];
rz(5.75419011910493) q[2];
sx q[2];
rz(6.44065163134738) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.46865141391754) q[3];
sx q[3];
rz(2.39904377062852) q[3];
sx q[3];
rz(9.31915531157657) q[3];
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
