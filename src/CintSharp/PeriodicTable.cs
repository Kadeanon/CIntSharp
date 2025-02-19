
/* 
 * Copyright (C) 2018  Kazuya Ujihara <ujihara.kazuya@gmail.com>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT Any WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */


using System;
using System.Collections.Generic;
using System.Collections.Immutable;

namespace CintSharp
{
    public record struct ChemicalElement(int AtomicNumber, string Symbol, string Name = "")
    {
        public override readonly string ToString()
        {
            return $"{AtomicNumber}:{Symbol ?? "null"}";
        }

        public static ChemicalElement Of(int number)
        {
            if (!(number >= 0 && number <= 118))
                ArgumentException.ThrowIfNullOrEmpty(null, nameof(number));
            return PeriodicTable.ElementDict[number];
        }

        public static ChemicalElement OfSymbol(string symbol)
        {
            return PeriodicTable.ElementDict.FirstOrDefault(
                element => element.Symbol == symbol, 
                PeriodicTable.Unknown);
        }

        public static ChemicalElement OfName(string name)
        {
            return PeriodicTable.ElementDict.FirstOrDefault(
                element => element.Name == name, 
                PeriodicTable.Unknown);
        }
    }

    public static class PeriodicTable
    {

        public static ImmutableArray<ChemicalElement> ElementDict { get; }

        static PeriodicTable()
        {
            ElementDict =
            [
                new ChemicalElement(0, "R", "Unknown"),
                new ChemicalElement(1, "H", "Hydrogen"),
                new ChemicalElement(2, "He", "Helium"),
                new ChemicalElement(3, "Li", "Lithium"),
                new ChemicalElement(4, "Be", "Beryllium"),
                new ChemicalElement(5, "B", "Boron"),
                new ChemicalElement(6, "C", "Carbon"),
                new ChemicalElement(7, "N", "Nitrogen"),
                new ChemicalElement(8, "O", "Oxygen"),
                new ChemicalElement(9, "F", "Fluorine"),
                new ChemicalElement(10, "Ne", "Neon"),
                new ChemicalElement(11, "Na", "Sodium"),
                new ChemicalElement(12, "Mg", "Magnesium"),
                new ChemicalElement(13, "Al", "Aluminium"),
                new ChemicalElement(14, "Si", "Silicon"),
                new ChemicalElement(15, "P", "Phosphorus"),
                new ChemicalElement(16, "S", "Sulfur"),
                new ChemicalElement(17, "Cl", "Chlorine"),
                new ChemicalElement(18, "Ar", "Argon"),
                new ChemicalElement(19, "K", "Potassium"),
                new ChemicalElement(20, "Ca", "Calcium"),
                new ChemicalElement(21, "Sc", "Scandium"),
                new ChemicalElement(22, "Ti", "Titanium"),
                new ChemicalElement(23, "V", "Vanadium"),
                new ChemicalElement(24, "Cr", "Chromium"),
                new ChemicalElement(25, "Mn", "Manganese"),
                new ChemicalElement(26, "Fe", "Iron"),
                new ChemicalElement(27, "Co", "Cobalt"),
                new ChemicalElement(28, "Ni", "Nickel"),
                new ChemicalElement(29, "Cu", "Copper"),
                new ChemicalElement(30, "Zn", "Zinc"),
                new ChemicalElement(31, "Ga", "Gallium"),
                new ChemicalElement(32, "Ge", "Germanium"),
                new ChemicalElement(33, "As", "Arsenic"),
                new ChemicalElement(34, "Se", "Selenium"),
                new ChemicalElement(35, "Br", "Bromine"),
                new ChemicalElement(36, "Kr", "Krypton"),
                new ChemicalElement(37, "Rb", "Rubidium"),
                new ChemicalElement(38, "Sr", "Strontium"),
                new ChemicalElement(39, "Y", "Yttrium"),
                new ChemicalElement(40, "Zr", "Zirconium"),
                new ChemicalElement(41, "Nb", "Niobium"),
                new ChemicalElement(42, "Mo", "Molybdenum"),
                new ChemicalElement(43, "Tc", "Technetium"),
                new ChemicalElement(44, "Ru", "Ruthenium"),
                new ChemicalElement(45, "Rh", "Rhodium"),
                new ChemicalElement(46, "Pd", "Palladium"),
                new ChemicalElement(47, "Ag", "Silver"),
                new ChemicalElement(48, "Cd", "Cadmium"),
                new ChemicalElement(49, "In", "Indium"),
                new ChemicalElement(50, "Sn", "Tin"),
                new ChemicalElement(51, "Sb", "Antimony"),
                new ChemicalElement(52, "Te", "Tellurium"),
                new ChemicalElement(53, "I", "Iodine"),
                new ChemicalElement(54, "Xe", "Xenon"),
                new ChemicalElement(55, "Cs", "Caesium"),
                new ChemicalElement(56, "Ba", "Barium"),
                new ChemicalElement(57, "La", "Lanthanum"),
                new ChemicalElement(58, "Ce", "Cerium"),
                new ChemicalElement(59, "Pr", "Praseodymium"),
                new ChemicalElement(60, "Nd", "Neodymium"),
                new ChemicalElement(61, "Pm", "Promethium"),
                new ChemicalElement(62, "Sm", "Samarium"),
                new ChemicalElement(63, "Eu", "Europium"),
                new ChemicalElement(64, "Gd", "Gadolinium"),
                new ChemicalElement(65, "Tb", "Terbium"),
                new ChemicalElement(66, "Dy", "Dysprosium"),
                new ChemicalElement(67, "Ho", "Holmium"),
                new ChemicalElement(68, "Er", "Erbium"),
                new ChemicalElement(69, "Tm", "Thulium"),
                new ChemicalElement(70, "Yb", "Ytterbium"),
                new ChemicalElement(71, "Lu", "Lutetium"),
                new ChemicalElement(72, "Hf", "Hafnium"),
                new ChemicalElement(73, "Ta", "Tantalum"),
                new ChemicalElement(74, "W", "Tungsten"),
                new ChemicalElement(75, "Re", "Rhenium"),
                new ChemicalElement(76, "Os", "Osmium"),
                new ChemicalElement(77, "Ir", "Iridium"),
                new ChemicalElement(78, "Pt", "Platinum"),
                new ChemicalElement(79, "Au", "Gold"),
                new ChemicalElement(80, "Hg", "Mercury"),
                new ChemicalElement(81, "Tl", "Thallium"),
                new ChemicalElement(82, "Pb", "Lead"),
                new ChemicalElement(83, "Bi", "Bismuth"),
                new ChemicalElement(84, "Po", "Polonium"),
                new ChemicalElement(85, "At", "Astatine"),
                new ChemicalElement(86, "Rn", "Radon"),
                new ChemicalElement(87, "Fr", "Francium"),
                new ChemicalElement(88, "Ra", "Radium"),
                new ChemicalElement(89, "Ac", "Actinium"),
                new ChemicalElement(90, "Th", "Thorium"),
                new ChemicalElement(91, "Pa", "Protactinium"),
                new ChemicalElement(92, "U", "Uranium"),
                new ChemicalElement(93, "Np", "Neptunium"),
                new ChemicalElement(94, "Pu", "Plutonium"),
                new ChemicalElement(95, "Am", "Americium"),
                new ChemicalElement(96, "Cm", "Curium"),
                new ChemicalElement(97, "Bk", "Berkelium"),
                new ChemicalElement(98, "Cf", "Californium"),
                new ChemicalElement(99, "Es", "Einsteinium"),
                new ChemicalElement(100, "Fm", "Fermium"),
                new ChemicalElement(101, "Md", "Mendelevium"),
                new ChemicalElement(102, "No", "Nobelium"),
                new ChemicalElement(103, "Lr", "Lawrencium"),
                new ChemicalElement(104, "Rf", "Rutherfordium"),
                new ChemicalElement(105, "Db", "Dubnium"),
                new ChemicalElement(106, "Sg", "Seaborgium"),
                new ChemicalElement(107, "Bh", "Bohrium"),
                new ChemicalElement(108, "Hs", "Hassium"),
                new ChemicalElement(109, "Mt", "Meitnerium"),
                new ChemicalElement(110, "Ds", "Darmstadtium"),
                new ChemicalElement(111, "Rg", "Roentgenium"),
                new ChemicalElement(112, "Cn", "Copernicium"),
                new ChemicalElement(113, "Nh", "Nihonium"),
                new ChemicalElement(114, "Fl", "Flerovium"),
                new ChemicalElement(115, "Mc", "Moscovium"),
                new ChemicalElement(116, "Lv", "Livermorium"),
                new ChemicalElement(117, "Ts", "Tennessine"),
                new ChemicalElement(118, "Og", "Oganesson"),
            ];
        }

        public static ChemicalElement Unknown => ElementDict[0];
    }
}
