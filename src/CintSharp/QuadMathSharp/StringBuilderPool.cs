using Microsoft.Extensions.ObjectPool;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace QuadMathSharp
{
    internal class StringBuilderPool : DefaultObjectPool<StringBuilder>
    {
        public StringBuilderPool(StringBuilderPooledObjectPolicy policy) : base(policy)
        {
        }

        public static StringBuilderPool Shared { get; }
            = new StringBuilderPool(new() { InitialCapacity = 255 });

        public StringBuilder Get(int capacity) 
        {
            var sb = base.Get();
            sb.EnsureCapacity(capacity);
            return sb;
        }

        public static PooledStringBuilderScope ScopedBorrow(int capacity, out StringBuilder entry)
        {
            return new PooledStringBuilderScope(capacity).Get(out entry);
        }
    }

    internal ref struct PooledStringBuilderScope : IDisposable 
    {
        private readonly ObjectPool<StringBuilder> _pool;
        private readonly StringBuilder _entry;
        private bool _disposed;
        public PooledStringBuilderScope(ObjectPool<StringBuilder>? pool = null)
        {
            _pool = pool??StringBuilderPool.Shared;
            _entry = _pool.Get();
            _disposed = false;
        }

        public PooledStringBuilderScope(int capacity, ObjectPool<StringBuilder>? pool = null)
        {
            _pool = pool ?? StringBuilderPool.Shared;
            _entry = _pool.Get();
            _entry.EnsureCapacity(capacity);
            _disposed = false;
        }

        public readonly PooledStringBuilderScope Get(out StringBuilder entry)
        {
            entry = _entry;
            return this;
        }

        public readonly StringBuilder? Entry => _entry;
        public void Dispose()
        {
            if(!_disposed)
            {
                _pool.Return(_entry);
                _disposed = true;
            }
        }
    }
}
